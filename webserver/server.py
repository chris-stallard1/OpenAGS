from quart import Quart, request, redirect, url_for, render_template, send_file, websocket, make_response
import asyncio
from hypercorn.middleware import HTTPToHTTPSRedirectMiddleware
from werkzeug.utils import secure_filename
import aiofiles
import aiofiles.os

import os
import uuid
import json
import time
import sys
import functools
import shutil
from copy import deepcopy

sys.path.append(os.path.join("..","backend"))

from backend import ActivationAnalysis
from constants import som
from evaluators import MassSensEval

loop = None

app = Quart(__name__)

activeProjects = dict()

@app.before_serving
async def startup():
    global loop
    loop = asyncio.get_event_loop() #globally store the event loop, so we can use it later

@app.route("/icons/<icon_name>")
async def get_icon(icon_name):
    """Gets the icon with name icon_name"""
    return await send_file(os.path.join(os.getcwd(), "icons", icon_name))

#Static Pages
@app.route("/")
async def homepage():
    async with aiofiles.open(os.path.join(os.getcwd(),'homepage.html'), mode='r') as f:
        contents = await f.read()
    return contents

@app.route("/error")
async def errorpage():
    async with aiofiles.open(os.path.join(os.getcwd(),'error.html'), mode='r') as f:
        contents = await f.read()
    return contents

#Dynamic but not Jinja Templates

@app.route("/create", methods=["GET","POST"])
async def create():
    """Handles GET requests (static page) and POST requests (file uploads/form submissions) for the create project page"""
    if request.method == "GET":
        async with aiofiles.open(os.path.join(os.getcwd(),'create.html'), mode='r') as f:
            contents = await f.read()
        return contents
    else:
        #POST Method, handle file uploads
        upload_path = os.path.join(os.getcwd(), "uploads")
        projectID = str(uuid.uuid4())
        #TODO: aiofiles.os.mkdir not working for some reason. Debug this
        
        os.mkdir(os.path.join(upload_path,projectID))
        os.mkdir(os.path.join(os.getcwd(), "results", projectID))

        uploaded_files = await request.files  
        form = await request.form      
        files_list = uploaded_files.getlist("file")
        try:
            standardsFile = uploaded_files.get("standardsFile")
            standardsFilename = os.path.join(upload_path,projectID,secure_filename(standardsFile.filename))
            await standardsFile.save(standardsFilename)
        except:
            #use the default file
            standardsFilename = os.path.join(os.getcwd(),"AllSensitivity.csv")
    
        filenamesList = []
        for f in files_list:
            filenamesList.append(os.path.join(upload_path,projectID,secure_filename(f.filename)))
            p = os.path.join(upload_path,projectID,secure_filename(f.filename))
            await f.save(p)
        
        #check analysis type
        delayed = form.get("analysisType") == "delayed"

        #initialize state file
        async with aiofiles.open(os.path.join(os.getcwd(), "uploads", projectID, "state.json"), mode="w") as f:
            await f.seek(0)
            await f.write(json.dumps({
                "title" : form["title"],
                "files" : filenamesList,
                "standardsFilename" : standardsFilename,
                "ROIsFitted" : False,
                "ROIs" : [],
                "resultsGenerated" : False,
                "delayed" : delayed,
                "NAATimes" : [[] for i in range(len(filenamesList))]
            }))
        return json.dumps({"id" : projectID})

@app.route("/results/<projectID>/<filename>")
async def serve_result(projectID, filename):
    """Serves a file from the results folder, containing the output of the program. 
    
    Basically a wrapper around the ActivationAnalysis.write_results_file() function that also serves the file if it has already been generated.
    """
    try:
        #if the file has already been generated, just serve it
        async with aiofiles.open(os.path.join(os.getcwd(),"results",projectID,filename), mode="rb") as f:
            c = await f.read()
            return c, 200, {'Content-Disposition' : 'attachment; filename="'+filename+'"'}
    except:
        #otherwise, generate it
        global activeProjects
        if projectID in activeProjects.keys(): #if its loaded in memory
            analysisObject = activeProjects[projectID]["analysisObject"]
        else: #load from state
            async with aiofiles.open(os.path.join(os.getcwd(),"uploads", projectID, "state.json"), mode="r") as f:
                contents = await f.read()
                currentProject = json.loads(contents)
                analysisObject = ActivationAnalysis()
                await loop.run_in_executor(None, analysisObject.load_from_dict, currentProject)
        analysisObject.write_results_file(projectID, filename)
        async with aiofiles.open(os.path.join(os.getcwd(),"results",projectID, filename), mode="rb") as f:
            c = await f.read()
            return c, 200, {'Content-Disposition' : 'attachment; filename="'+filename+'"'}

#Jinja Template Pages

@app.route("/projects/<projectID>/<action>")
async def project(projectID, action):
    global activeProjects
    analysisObject = None
    if projectID in activeProjects.keys(): #if the project is loaded
        analysisObject = activeProjects[projectID]["analysisObject"]
    else: #load the project from the state.json file
        async with aiofiles.open(os.path.join(os.getcwd(),"uploads",projectID,"state.json"), mode="r") as f:
            contents = await f.read()
            currentProject = json.loads(contents)
            analysisObject = ActivationAnalysis()
            await loop.run_in_executor(None, analysisObject.load_from_dict, currentProject)
            activeProjects[projectID] = {
                "analysisObject" : analysisObject,
                "webSockets" : [],
                "numUsers" : 0,
                "saveAction" : None
            }
    
    if action == "edit":
        if not analysisObject.ROIsFitted:
            analysisObject.get_fitted_ROIs()
        return await(render_template("project.html", analysisObject=analysisObject, projectID=projectID, pathSplit=os.path.split))
    elif action == "view":
        return await(render_template("view.html", analysisObject=analysisObject, projectID=projectID, som=som, pathSplit=os.path.split))
    elif action == "results":
        return await(render_template("results.html", analysisObject=analysisObject, projectID=projectID, pathSplit=os.path.split))
    elif action == "delete":
        p1 = os.path.join(os.getcwd(), "results", projectID)
        if os.path.exists(p1):
            shutil.rmtree(p1)
        p2 = os.path.join(os.getcwd(), "uploads", projectID)
        if os.path.exists(p2):
            shutil.rmtree(p2)
        activeProjects[projectID]["saveAction"].cancel()
        del activeProjects[projectID]
        return redirect("/create")
    else: 
        return ""

#WebSocket handler Utility Functions

async def saveProject(projectID):
    """Save the project after 60 seconds of a user being off of the page"""
    await asyncio.sleep(60)
    print("saving "+projectID)
    global activeProjects
    if projectID not in activeProjects.keys():
        return None
    async with aiofiles.open(os.path.join(os.getcwd(),"uploads",projectID,"state.json"), mode="w") as f:
        await f.seek(0)
        await f.write(json.dumps(activeProjects[projectID]["analysisObject"].export_to_dict()))
    
    del activeProjects[projectID]

async def saveProjectNow(projectID):
    """Quicksave function"""
    print("saving "+projectID)
    global activeProjects
    if projectID not in activeProjects.keys():
        return None
    async with aiofiles.open(os.path.join(os.getcwd(),"uploads",projectID,"state.json"), mode="w") as f:
        await f.seek(0)
        await f.write(json.dumps(activeProjects[projectID]["analysisObject"].export_to_dict()))

#WebSocket Handler function

@app.websocket("/projects/<projectID>/ws")
async def ws(projectID):
    """The WebSocket Handler Function
    
    Producer/Consumer Functions explained at: https://medium.com/@pgjones/websockets-in-quart-f2067788d1ee under "Consuming and Producing Tasks"
    Use of asyncio.Queue for broadcast taken from: https://pgjones.gitlab.io/quart/tutorials/websocket_tutorial.html under "Section 5: Broadcasting"
    """

    async def producer(projectID):
        global activeProjects
        queue = asyncio.Queue()
        if projectID not in activeProjects.keys():
            #load project
            async with aiofiles.open(os.path.join(os.getcwd(),"uploads",projectID,"state.json"), mode="r") as f:
                contents = await f.read()
                currentProject = json.loads(contents)
                analysisObject = ActivationAnalysis()
                await loop.run_in_executor(None, analysisObject.load_from_dict, currentProject)
                activeProjects[projectID] = {
                    "analysisObject" : analysisObject,
                    "webSockets" : [],
                    "numUsers" : 0,
                    "saveAction" : None
                }
        activeProjects[projectID]["webSockets"].append(queue)
        activeProjects[projectID]["numUsers"] += 1
        if activeProjects[projectID]["saveAction"] != None:#if we were going to save and remove the project, don't
            activeProjects[projectID]["saveAction"].cancel()
        activeProjects[projectID]["saveAction"] = None
        print("save action cancelled")
        while True:
            try:#listen on created queue, this works for broadcasting
                data = await queue.get()
                await websocket.send(data)
            except asyncio.CancelledError:
                pass

    async def consumer(projectID):
        global activeProjects
        while True:
            data = await websocket.receive()
            dataDict = json.loads(data)
            if dataDict["type"] == "ROIUpdate":
                #Refit an ROI
                analysisObject = activeProjects[projectID]["analysisObject"]
                ROI = analysisObject.ROIs[dataDict["index"]]
                ROI.fitted = False
                if "newRange" in dataDict.keys():
                    analysisObject.set_ROI_range(dataDict["index"],dataDict["newRange"])
                #remove the peaks that were removed by the user
                if "existingPeaks" in dataDict.keys():
                    peaks = ROI.get_peaks()
                    newPeaks = []
                    for peak in peaks:
                        if peak.to_string() in dataDict["existingPeaks"]:
                            newPeaks.append(peak)
                    ROI.set_peaks(newPeaks)
                #add peaks added by the user
                if "newPeaks" in dataDict.keys():
                    for peak in dataDict["newPeaks"]:
                        peakType = peak[0]
                        peakParams = peak[1:]
                        ROI.peaks.append(som["peaks"][peakType](*peakParams))
                #change the background if necessary
                if "background" in dataDict.keys():
                    bg = dataDict["background"]
                    bgType = bg[0]
                    bgParams = bg[1:]
                    ROI.set_background(som["backgrounds"][bgType](*bgParams))
                #call the fitter
                ROI.fit()
                #prepare the output
                outputObj = {
                    "type" : "ROIUpdate",
                    "index" : dataDict["index"],
                    "fitted" : ROI.fitted,
                    "xdata" : ROI.get_energies(),
                    "ydata" : ROI.get_cps(),
                    "peakStrings" : [p.to_string() for p in ROI.get_peaks()],
                    "bgString" : ROI.get_background().to_string(),
                    "ROIRange" : ROI.get_range()
                    }
                if ROI.fitted:
                    curve = ROI.get_fitted_curve()
                    outputObj["curveX"] = curve[0]
                    outputObj["curveY"] = curve[1]
                    outputObj["peakX"] = ROI.get_peak_ctrs()
                    outputObj["backgroundY"] = curve[2]

                data = json.dumps(outputObj)   

                #Broadcast
                for queue in activeProjects[projectID]["webSockets"]:
                    await queue.put(data)
            
            elif dataDict["type"] == "entryReprRequest":
                analysisObject = activeProjects[projectID]["analysisObject"]
                entryParams = []
                try:#test if inputs are floats
                    entryParams = [float(x) for x in dataDict["entryParams"]]
                except:
                    await websocket.send(json.dumps({"type" : "error", "text":"Please enter only numbers for peak paramaters."}))
                    continue
                try:#test if inputs are within  ROI bounds, and get results if so
                    stringRepr, params = analysisObject.get_entry_repr(dataDict["class"],dataDict["name"],dataDict["ROIIndex"],entryParams)
                except:
                    await websocket.send(json.dumps({"type" : "error", "text":"Your peak is outside the ROI bounds."}))
                    continue #must continue instead of breaking because breaking kills the websocket
                #prepare output
                outputObj = {
                    "type" : "entryReprResponse",
                    "class" : dataDict["class"],
                    "index" : dataDict["ROIIndex"],
                    "name" : dataDict["name"],
                    "params" : params,
                    "output" : stringRepr
                }
                #Broadcast
                for queue in activeProjects[projectID]["webSockets"]:
                    await queue.put(json.dumps(outputObj))
            elif dataDict["type"] == "matchUpdate":
                #Update peak matches by echoing to all users
                for queue in activeProjects[projectID]["webSockets"]:
                        await queue.put(data)

            elif dataDict["type"] == "titleUpdate":
                #Update the title
                analysisObject = activeProjects[projectID]["analysisObject"]
                if dataDict["newTitle"] != analysisObject.get_title(): #make sure it changed
                    analysisObject.set_title(dataDict["newTitle"])
                    #Broadcast
                    for queue in activeProjects[projectID]["webSockets"]:
                        await queue.put(data)
            elif dataDict["type"] == "NAATimeUpdate":
                #Update NAA Times (irradiation, wait time)
                analysisObject = activeProjects[projectID]["analysisObject"]
                analysisObject.fileData[dataDict["fileIndex"]]["NAATimes"] = dataDict["times"]
                #Broadcast
                for queue in activeProjects[projectID]["webSockets"]:
                    await queue.put(data)
            elif dataDict["type"] == "isotopeUpdate":
                #Update isotopes 
                analysisObject = activeProjects[projectID]["analysisObject"]
                analysisObject.update_ROIs(dataDict["addedIsotopes"], dataDict["removedIsotopes"])
                outputObj = {
                    "type" : "isotopeUpdateResponse",
                    "currentIsotopes" : analysisObject.get_isotopes(),
                    "ROIRanges" : [r.get_formatted_range() for r in analysisObject.ROIs],
                    "ROIIndicies" : [r.get_indicies() for r in analysisObject.ROIs],
                    "ROIIsotopes" : [', '.join(r.get_isotopes()) for r in analysisObject.ROIs]
                }
                outputData = json.dumps(outputObj)
                #Broadcast
                for queue in activeProjects[projectID]["webSockets"]:
                    await queue.put(outputData)
            elif dataDict["type"] == "peakMatchSubmission":
                #Generate Results
                analysisObject = activeProjects[projectID]["analysisObject"]
                for i in range(len(dataDict["pairs"])):
                    analysisObject.ROIs[i].set_original_peak_pairs(dataDict["pairs"][i])
                #TODO: Maybe allow users to customize this, as that is kinda the point of evaluators. 
                analysisObject.run_evaluators([MassSensEval], [[]])
                
                await saveProjectNow(projectID)
                
                for f in os.listdir(os.path.join(os.getcwd(), "results", projectID)):
                    os.remove(os.path.join(os.getcwd(), "results", projectID, f))

                #Broadcast simple message which will redirect people to the results screen
                outputData = json.dumps({"type" : "resultsGenerated"})
                for queue in activeProjects[projectID]["webSockets"]:
                    await queue.put(outputData)

            elif dataDict["type"] == "userPrefsUpdate":
                #Update User Preferences (settings)
                analysisObject = activeProjects[projectID]["analysisObject"]
                analysisObject.set_user_prefs(dataDict["newPrefs"])
                for queue in activeProjects[projectID]["webSockets"]:
                    await queue.put(data)
    consumer_task = asyncio.ensure_future(consumer(projectID))
    producer_task = asyncio.ensure_future(producer(projectID))
    try:
        await asyncio.gather(consumer_task, producer_task)
    finally:
        activeProjects[projectID]["numUsers"] -= 1
        if activeProjects[projectID]["numUsers"] <= 0: # save project in 1 min if no one reconnects
            print("started save action")
            activeProjects[projectID]["saveAction"] = asyncio.create_task(saveProject(projectID))
        consumer_task.cancel()
        producer_task.cancel()
@app.after_serving
async def export_to_db():
    global activeProjects
    #Save all the projects
    for projectID in activeProjects.keys():
        a = activeProjects[projectID]["saveAction"]
        if a != None:
            a.cancel()
        await saveProjectNow(projectID)

redirectedApp = HTTPToHTTPSRedirectMiddleware(app, host="openags.cstallar.net")