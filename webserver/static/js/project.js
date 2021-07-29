//unintended but incredibly smart solution that uses ws:// if we are on http, but wss:// if we are on https
var wsUrl = window.location.href.toString().replace("http","ws").replace("edit","ws");
var ws = new WebSocket(wsUrl);
ws.onmessage = function (event) {
    var data = JSON.parse(event.data);
    switch(data.type){
        case "titleUpdate":
            //just update the title
            document.getElementById("cdTitle").value = data.newTitle;
            break;
        case "entryReprResponse":
            //add an entry to the current fit section
            if(data.class == "peaks"){
                var peaksList = document.getElementById("fittedPeaksList-"+data.index.toString());
                var output = data.output;
                var peakData = [data.name].concat(data.params);
                window.newPeaks[data.index][output] = peakData;
                peaksList.innerHTML = peaksList.innerHTML + "<li class='list-group-item compact' id=\""+output+"\"><p class = 'peak-label'>"+output+"</p><img class='rmv-btn' src='/icons/file-x.svg' onclick='remove_peak_from_list(\""+ output +"\")'/></li>";
            }
            else if(data.class == "backgrounds"){
                var bgList = document.getElementById("fittedBgList-"+data.index.toString());
                var output = data.output;
                var bgData = [data.name].concat(data.params);
                window.newBackgrounds[data.index] = bgData;
                bgList.innerHTML = DOMPurify.sanitize("<li class='list-group-item compact' id=\""+output+"\"><p class = 'bg-label'>"+output+"</p>");
            }
            break;
        case "ROIUpdate":
            //update absolutely everything, new fit
            var i = data.index;

            //reset csome stuff
            window.newPeaks[i] = {};
            window.newBackgrounds[i] = null;
            window.originalEnergyBounds[i] = data.ROIRange;
            resetEnergyBounds(i);

            dataTrace = {
                x : data.xdata,
                y : data.ydata,
                mode : "markers",
                name : "Data"
            }
            //if we have fitted a previously unfitted region, remove it from unfitetd regions list
            if(data.fitted && window.unfittedRegions.indexOf(data.index) !== -1){
                window.unfittedRegions.splice(window.unfittedRegions.indexOf(data.index), 1);
            }

            //if we have messed up a previously fitted region, add it to unfitted list and show an error
            if(!data.fitted){
                if(window.unfittedRegions.indexOf(data.index) === -1){
                    window.unfittedRegions.push(data.index);
                }
                showErrorMessage("Could not find a fit. Try reducing the number of peaks in the fit. If this doesn't work, you can submit the other ROIs and this one will be ignored.");
                break;
            }

            //set new background
            var bgList = document.getElementById("fittedBgList-"+data.index.toString());
            bgList.innerHTML = DOMPurify.sanitize("<li class='list-group-item compact' id=\""+ data.bgString +"\"><p class = 'bg-label'>"+ data.bgString +"</p>");
            
            //set new peaks
            var peaksList = document.getElementById("fittedPeaksList-"+data.index.toString());
            peaksList.innerHTML = "";
            for(j=0;j<data.peakStrings.length;j++){
                var output = data.peakStrings[j];
                peaksList.innerHTML = peaksList.innerHTML + "<li class='list-group-item compact' id=\""+output+"\"><p class = 'peak-label'>"+output+"</p><img class='rmv-btn' src='/icons/file-x.svg' onclick='remove_peak_from_list(\""+ output +"\")'/></li>";
            }

            //set new matches
            var knownPeaks = document.getElementById("userPeakMatches-"+i.toString()).getElementsByTagName("p");
            var matchSelects = document.getElementById("userPeakMatches-"+i.toString()).getElementsByClassName("peak-match-select");
            for(j=0;j<matchSelects.length;j++){
                var selectHTML = '';
                var knownCtr = parseFloat(knownPeaks[j].id);
                var minSep = 99;
                for(k=0;k<data.peakX.length;k++){
                    if(Math.abs(data.peakX[k] - knownCtr) < minSep){
                        minSep = Math.abs(data.peakX[k] - knownCtr);
                    }
                }
                for(k=0;k<data.peakStrings.length;k++){
                    if(Math.abs(data.peakX[k] - knownCtr) === minSep){
                        selectHTML += "<option selected value='"+data.peakX[k].toString()+"'>"+data.peakStrings[k]+"</option>"
                    }
                    else{
                        selectHTML += "<option value='"+data.peakX[k].toString()+"'>"+data.peakStrings[k]+"</option>"
                    }
                }
                matchSelects[j].innerHTML = DOMPurify.sanitize(selectHTML);
            }

            var plot = document.getElementById("ROI-"+i.toString());

            //add other traces if the data is fitted, and update the plot
            if(data.fitted){
                var fitTrace = {
                    x : data.curveX,
                    y : data.curveY,
                    mode : "lines",
                    name : "Fit"
                }
                var peakY = data.peakX.map(x => data.curveY[Math.floor((x - data.curveX[0])*100)]);
                var peakTrace = {
                    x : data.peakX,
                    y : peakY,
                    mode : "markers",
                    name : "Peaks" 
                }
                var bgTrace = {
                    x : data.curveX,
                    y : data.backgroundY,
                    mode : "lines",
                    name : "Background"
                }
                var allTraces = [bgTrace, fitTrace, peakTrace, dataTrace];
                while(plot.data.length > 0){
                    Plotly.deleteTraces(plot, 0);
                }
                Plotly.addTraces(plot, allTraces);
            }
            else{
                var allTraces = [dataTrace];
                while(plot.data.length > 0){
                    Plotly.deleteTraces(plot, 0);
                }
                Plotly.addTraces(plot, allTraces);
            }
            break;
        case "matchUpdate":
            //update one of the matches
            document.getElementById("userPeakMatches-"+data.ROIIndex.toString()).getElementsByTagName('select')[data.matchIndex.toString()].value = data.newValue.toString();
            break;
        case "isotopeUpdateResponse":
            //if isotopes are updated, the page must be refreshed
            ws.close();       
            showModal("refreshPageModal");
            break;
        case "resultsGenerated":
            //if results are generated, you usualyl want to view them
            ws.close();
            showModal("redirectModal");
            break;
        case "error":
            //if there's a backend error, show it in the frontend
            showErrorMessage(data.text);
            break;
    }

};

/**
 * Sets the peak entry options below the "Select Peak" box in the ith ROI window. uses window.entryFields to see field names and count.
 * @param {Number} i 
 */
function updatePeakEntry(i){
    var peakType = document.getElementById("peakSelect-"+i.toString()).value;
    var peaksList = document.getElementById("userPeakEntry-"+i.toString());
    var newInnerHTML = "";
    if(peakType !== "Select Peak Type"){
        window.entryFields["peaks"][peakType].forEach((field) => {newInnerHTML = newInnerHTML + "<li class='list-group-item'><p style='float:left;'>"+field+"</p><input class='form-control w-50 peak-entry' style='float:right;'/></li>"});
    }
    peaksList.innerHTML = DOMPurify.sanitize(newInnerHTML);
};

/**
 * Transmits an added peak through the WebSocket for ROI #i
 * @param {Number} i 
 */
function add_peak_to_list(i){
    var peakType = document.getElementById("peakSelect-"+i.toString()).value;
    var peaksList = document.getElementById("userPeakEntry-"+i.toString());
    var peakProps = peaksList.getElementsByClassName("peak-entry");
    var propValues = [];
    for(j=0;j<peakProps.length;j++){
        propValues[j] = peakProps[j].value;
        peakProps[j].value = "";
    }
    var wsSendObj = {
        "type" : "entryReprRequest",
        "ROIIndex" : i,
        "class" : "peaks",
        "name" : peakType,
        "entryParams" : propValues
    }
    ws.send(JSON.stringify(wsSendObj));
};

//TODO: Potentially unnecessary, look at removing
function remove_peak_from_list(peakID){
    document.getElementById(peakID).remove();
}

/**
 * Transmits a user-selected background type to the ebsocket for ROI #i
 * @param {Number} i 
 */
function editBackground(i){
    var bgType = document.getElementById("backgroundSelect-"+i.toString()).value;
    var wsSendObj = {
        "type" : "entryReprRequest",
        "ROIIndex" : i,
        "class" : "backgrounds",
        "name" : bgType,
        "entryParams" : []
    }
    ws.send(JSON.stringify(wsSendObj));
}

/**
 * Resets the ROI range for ROI #i in case a user has entered an invalid value
 * @param {Number} i 
 */
function resetEnergyBounds(i){
    var entryList = document.getElementById("editRangeList-"+i.toString()).getElementsByTagName("input");
    entryList[0].value = window.originalEnergyBounds[i][0].toString();
    entryList[1].value = window.originalEnergyBounds[i][1].toString();
}
/**
 * Sends the WebSocket request to renaalyze ROI #i
 * @param {Number} i 
 */

function reanalyze(i){
    var peaksList = document.getElementById("fittedPeaksList-"+i.toString());
    var peaks = peaksList.getElementsByClassName("peak-label");
    
    var existingPeaksToKeep = [];
    var newPeaksToAdd = [];
    for(j=0;j<peaks.length;j++){
        peak = peaks[j].innerText;
        if(window.newPeaks[i].hasOwnProperty(peak)){//check if this peak was already added last time
            newPeaksToAdd.push(window.newPeaks[i][peak]);
        }
        else{
            existingPeaksToKeep.push(peak);
        }
    }
    var bgToAdd = window.newBackgrounds[i];
    
    outputObject = {"type" : "ROIUpdate","index":i};

    //add new range to the object, if applicable 
    var entryList = document.getElementById("editRangeList-"+i.toString()).getElementsByTagName("input");
    if(entryList[0].value != window.originalEnergyBounds[i][0] || entryList[1].value != window.originalEnergyBounds[i][1]){
        var minEnergy = parseFloat(entryList[0].value);
        var maxEnergy = parseFloat(entryList[1].value);
        if(isNaN(minEnergy) || isNaN(maxEnergy)){
            return showErrorMessage("Please enter decimal numbers for the range.")
        }
        outputObject["newRange"] = [minEnergy, maxEnergy];
    }

    //add peaks to keep, add, plus the background to the request
    if(existingPeaksToKeep !== []){
        outputObject["existingPeaks"] = existingPeaksToKeep;
    }
    if(newPeaksToAdd !== []){
        outputObject["newPeaks"] = newPeaksToAdd;
    }
    if(bgToAdd !== null){
        outputObject["background"] = bgToAdd;
    }
    console.log(JSON.stringify(outputObject));
    ws.send(JSON.stringify(outputObject));
}

function updateTitle(){
    var newTitle = document.getElementById("cdTitle").value;
    ws.send('{"type":"titleUpdate","newTitle":"'+newTitle+'"}');
}

/**
 * Submits the user's matches of known peaks and peaks in data, which will then bring up the results screen.
 */
function submitMatches(){
    var peakPairs = [];
    for(var i=0;i<numberPages;i++){
        var numKnownPeaks = document.getElementById("userPeakMatches-"+i.toString()).childElementCount;
        var knownPeakLabels = document.getElementById("userPeakMatches-"+i.toString()).getElementsByTagName("p");
        var peakMatches = document.getElementById("userPeakMatches-"+i.toString()).getElementsByTagName("select");
        var pairsInRegion = []
        for(var j=0;j<numKnownPeaks;j++){
            var pair = [parseFloat(peakMatches[j].value), parseFloat(knownPeakLabels[j].id)];
            pairsInRegion.push(pair);
        }
        peakPairs.push(pairsInRegion);
    }
    var outputObj = {
        "type" : "peakMatchSubmission",
        "pairs" : peakPairs
    }
    ws.send(JSON.stringify(outputObj));
}

/**
 * Sends a WebSocket request saying that the user has updated the jth match of the ith ROI, and sends the new value.
 * @param {Number} i 
 * @param {Number} j 
 */
function sendMatchUpdate(i,j){
    var newValue = document.getElementById("userPeakMatches-"+i.toString()).getElementsByTagName('select')[j.toString()].value;
    var outputObj = {
        "ROIIndex" : i,
        "matchIndex" : j,
        "newValue" : newValue
    };
    ws.send(JSON.stringify(outputObj));
}

function confirmSubmission(){
    if(window.unfittedRegions.length > 0){
        document.getElementById("confirmWarningText").innerText = "Unfitted Regions: ("+window.unfittedRegions.map(x => x+1).toString()+")";
    }
    else{
        document.getElementById("confirmWarningText").innerText = "";
    }
    showModal('confirmSubmitModal');
}