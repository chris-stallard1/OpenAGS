//unintended but incredibly smart solution that uses ws:// if we are on http, but wss:// if we are on https
var wsUrl = window.location.href.toString().replace("http","ws").replace("view","ws");
var ws = new WebSocket(wsUrl);
ws.onmessage = function (event) {
    var data = JSON.parse(event.data);
    switch(data.type){
        case "titleUpdate":
            document.getElementById("cdTitle").value = data.newTitle;
            break;
        case "isotopeUpdateResponse":
            var currentIsotopes = data.currentIsotopes;
            var newInnerHTML = "";
            for(var i=0;i<currentIsotopes.length;i++){
                var isotope = currentIsotopes[i];
                var newElement = "<p class='iso-label'>" + isotope  + "</p>";
                var removeButton = '<img class="rmv-btn" src="/icons/file-x.svg" onclick="removeIsotope('+"'"+isotope+"'"+')">';
                newInnerHTML += "<li class='list-group-item w-50' id='"+isotope+"'>" + newElement + removeButton +"</li>";
            }
            document.getElementById("selectedIsotopes").innerHTML = DOMPurify.sanitize(newInnerHTML);
            var newSelectHTML = "";
            for(var i=0;i<data.ROIRanges.length;i++){
                newSelectHTML += "<option value='"+data.ROIRanges[i][0]+","+data.ROIRanges[i][1]+","+data.ROIIndicies[i][0].toString()+","+data.ROIIndicies[i][1].toString()+"'>"+data.ROIRanges[i][0]+"-"+data.ROIRanges[i][1]+" keV ("+data.ROIIsotopes[i]+")</option>"
            }
            
            document.getElementById("compRangeSelect").innerHTML = DOMPurify.sanitize("<option selected value='custom'>Custom Range</option>" + newSelectHTML);
            for(var i=0;i<filesList.length;i++){
                document.getElementById("zoomToRegion-"+i.toString()).innerHTML = DOMPurify.sanitize("<option selected value=''>Whole Spectrum</option>" + newSelectHTML);
            }

            break;
        case "NAATimeUpdate":
            NAATimes[data.fileIndex] = data.times;
            updateShownTimes();
            break;
        case "error":
            showErrorMessage(data.text);
            break;
        case "userPrefsUpdate":
            document.getElementById("prefROIWidth").value = data.roi_width.toString();
            document.getElementById("prefBoronROIWidth").value = data.B_roi_width.toString();
            document.getElementById("prefPeakType").value = data.peak_type;
            document.getElementById("prefBoronPeakType") = data.boron_peak_type;
            document.getElementById("prefBGType") = data.background_type;
            document.getElementById("overlapROICheck").checked = data.overlap_rois;
            break;
    }
};

/**
 * Checks to ensure that ROIs have been selected and (if applicable) that times have been added, then opens analysis page
 */
function startAnalysis(){
    if(document.getElementById("selectedIsotopes").childElementCount === 0){
        showModal("updateROIModal"); //make the user add ROIs if they haven't 
    }
    else{
        try {
            var sw = false;
            for(var i=0;i<NAATimes.length;i++){
                if(NAATimes[i].length === 0){
                    sw=true;
                    document.getElementById("timeFileSelect").value = filesList[i];
                    updateShownTimes();
                    showModal("timeEntryModal");//make the user add times if they haven't
                }
            }
            if(!sw){
                window.location.replace(window.location.href.replace("/view","/edit"));
            }
        } catch (error) {
            console.log(error);
            window.location.replace(window.location.href.replace("/view","/edit"));
        }
        
    }
}

/**
 * Send an update when the user presses the submit button in the ROI modal. This keeps the backend and other users syncd.
 */
function submitROIs(){
    var wsObj = {
        "type" : "isotopeUpdate",
        "addedIsotopes" : addedIsotopes,
        "removedIsotopes" : removedIsotopes
    }
    ws.send(JSON.stringify(wsObj));
}

/**
 * Send a WebSocket request to update the times for NAA
 */
function sendNAATimes(){
    var irrTime = parseFloat(document.getElementById("irrTimeInput").value);
    var waitTime = parseFloat(document.getElementById("waitTimeInput").value);
    if(isNaN(irrTime) || isNaN(waitTime)){
        return showErrorMessage("Please enter times as numbers, in minutes.");
    }
    var allTimes = [irrTime, waitTime];
    var wsObj = {
        "type" : "NAATimeUpdate",
        "fileIndex" : filesList.indexOf(document.getElementById("timeFileSelect").value),
        "times" : allTimes
    };
    ws.send(JSON.stringify(wsObj));
}

/**
 * Send a WebSocket request to update the analysis settings (user preferences)
 */
function sendPrefUpdates(){
    var ROIWidth = parseFloat(document.getElementById("prefROIWidth").value);
    var boronROIWidth = parseFloat(document.getElementById("prefBoronROIWidth").value);
    if(isNaN(ROIWidth) || isNaN(boronROIWidth) || ROIWidth < 0 || boronROIWidth < 0){
        return showErrorMessage("Please enter ROI widths as positive numbers.")
    }
    var peakType = document.getElementById("prefPeakType").value;
    var boronPeakType = document.getElementById("prefBoronPeakType").value;
    var bgType = document.getElementById("prefBGType").value;
    var overlapROIs = document.getElementById("overlapROICheck").checked;
    wsObj = {
        "type" : "userPrefsUpdate",
        "newPrefs" : {
            "roi_width" : ROIWidth,
            "B_roi_width" : boronROIWidth,
            "peak_type" : peakType,
            "boron_peak_type" : boronPeakType,
            "background_type" : bgType,
            "overlap_rois" : overlapROIs
        }
    }
    ws.send(JSON.stringify(wsObj));
}

/**
 * Updates the document title (sends WebSocket)
 */
function updateTitle(){
    var newTitle = document.getElementById("cdTitle").value;
    ws.send('{"type":"titleUpdate","newTitle":"'+newTitle+'"}');
}

/**
 * Update the graph in the compare modal, using all info the user has entered
 */
function updateCompareModal(){
    var filename1 = document.getElementById("file1Select").value;
    var filename2 = document.getElementById("file2Select").value;
    var file1Index = filesList.indexOf(filename1);
    var file2Index = filesList.indexOf(filename2);
    var rangeSelectValue = document.getElementById("compRangeSelect").value;
    if(rangeSelectValue == "custom"){
        var minEnergy = parseFloat(document.getElementById("lowerBoundInput").value);
        var maxEnergy = parseFloat(document.getElementById("upperBoundInput").value);
    }
    else{
        var range = rangeSelectValue.split(",")
        var minEnergy = parseFloat(range[0]);
        var maxEnergy = parseFloat(range[1]);
    }
    var plot1 = document.getElementById("file-"+file1Index);
    var plot2 = document.getElementById("file-"+file2Index);
    var xData1 = plot1.data[0].x.slice(findClosest(plot1.data[0].x,minEnergy), findClosest(plot1.data[0].x,maxEnergy));
    var xData2 = plot2.data[0].x.slice(findClosest(plot2.data[0].x,minEnergy), findClosest(plot2.data[0].x,maxEnergy));
    var yData1 = plot1.data[0].y.slice(findClosest(plot1.data[0].x,minEnergy), findClosest(plot1.data[0].x,maxEnergy));
    var yData2 = plot2.data[0].y.slice(findClosest(plot2.data[0].x,minEnergy), findClosest(plot2.data[0].x,maxEnergy));
    var overlayGraphs = document.getElementById("overlayCheckbox").checked;
    if(overlayGraphs){
        var data = [
            {
                x : xData1,
                y : yData1,
                mode : 'lines',
                name: filename1
            },
            {
                x: xData2,
                y: yData2,
                mode: 'lines',
                name: filename2
            }
        ];
        var layout = {
            title: "Comparison ("+minEnergy+"-"+maxEnergy+" keV)",
            showlegend: false,
            xaxis: {
                "title" : "Energy (keV)"
            },
            yaxis: {
                "title" : "Counts Per Second"
            }
        };
        Plotly.react(document.getElementById("compareSpectraPlot"),data,layout,universalPlotConfig);
    }
    else{
        var data = [
            {
                x : xData1,
                y : yData1,
                mode : 'lines',
                name: filename1,
                xaxis : "x",
                yaxis : "y"
            },
            {
                x: xData2,
                y: yData2,
                mode: 'lines',
                name: filename2,
                xaxis : "x2",
                yaxis : "y2"
            }
        ];
        var layout = {
            title: "Comparison ("+minEnergy+"-"+maxEnergy+" keV)",
            showlegend: false,
            xaxis: { 
              domain: [0,0.48] ,
              title: "Energy (keV)"
            },
            yaxis: { 
              domain: [0,1],
              title: "Counts Per Second"
            },
            xaxis2: {
              domain: [0.52, 1],
              anchor: "y2",
              title: "Energy (keV)"
            },
            yaxis2: {
              domain: [0, 1],
              anchor: "x2"
            },
            annotations: [
              {
                text: "File 1",
                showarrow: false,
                x: 0,
                xref: "x domain",
                y: 1.1,
                yref: "y domain"
              },
              {
                text: "File 2",
                showarrow: false,
                x: 0,
                xref: "x2 domain",
                y: 1.1,
                yref: "y2 domain"
              }
            ]
          };
          Plotly.react(document.getElementById("compareSpectraPlot"),data,layout,universalPlotConfig);
    }
}

/**
 * Zoom into a region of file #i
 * @param {Number} i 
 */
function zoomToRegion(i){
    var selectObject = document.getElementById("zoomToRegion-"+i.toString());
    if(selectObject.value === ""){
        var plot = document.getElementById("file-"+i.toString());
        var xdata = plot.data[0].x
        document.getElementById("minEnergyInput-"+i.toString()).value = xdata[0];
        document.getElementById("maxEnergyInput-"+i.toString()).value = xdata[xdata.length - 1];
        var newLayout = {
            xaxis : {
                title: plot.layout.xaxis.title,
                range: [xdata[0], xdata[xdata.length - 1]]
            },
            yaxis : {
                title: plot.layout.yaxis.title,
                type : plot.layout.yaxis.type,
                range : [0, Math.max(...plot.data[0].y)*1.1]
            }
        };
    }
    else{
        var values = selectObject.value.split(",");
        var plot = document.getElementById("file-"+i.toString());
        var dataRange = plot.data[0].y.slice(parseInt(values[2]), parseInt(values[3]));
        document.getElementById("minEnergyInput-"+i.toString()).value = values[0];
        document.getElementById("maxEnergyInput-"+i.toString()).value = values[1];
        var newLayout = {
            xaxis : {
                title: plot.layout.xaxis.title,
                range: [parseFloat(values[0]), parseFloat(values[1])]
            },
            yaxis : {
                title: plot.layout.yaxis.title,
                type : plot.layout.yaxis.type,
                range : [0, Math.max(...dataRange)*1.1]
            }
        };
    }
    Plotly.relayout(plot, newLayout);
}

/**
 * Update the range of file #i, based on manual entries in the range input boxes
 * @param {Number} i 
 */
function updateRange(i){
    var minEnergy = parseFloat(document.getElementById("minEnergyInput-"+i.toString()).value);
    var maxEnergy = parseFloat(document.getElementById("maxEnergyInput-"+i.toString()).value);
    var plot = document.getElementById("file-"+i.toString());
    if(isNaN(minEnergy) || isNaN(maxEnergy) || maxEnergy <= minEnergy){
        var range = plot.layout.xaxis.range;
        document.getElementById("minEnergyInput-"+i.toString()).value = range[0].toString();
        document.getElementById("maxEnergyInput-"+i.toString()).value = range[1].toString();
        return showErrorMessage("Please enter decimal numbers for the data range, with Max. Energy > Min. Energy.");
    }
    var lowerIndex = findClosest(plot.data[0].x, minEnergy);
    var upperIndex = findClosest(plot.data[0].x, maxEnergy);
    var dataRange = plot.data[0].y.slice(lowerIndex, upperIndex);
    var myLayout = {
        xaxis : {
            title: plot.layout.xaxis.title,
            range: [minEnergy, maxEnergy]
        },
        yaxis : {
            title: plot.layout.yaxis.title,
            type : plot.layout.yaxis.type,
            range : [0, Math.max(...dataRange)*1.1]
        }
    };
    Plotly.relayout(plot, myLayout);
}

addedIsotopes = [];
removedIsotopes = [];

/**
 * Add an isotope to the analysis
 * @param {String} isotope 
 */
function addIsotope(isotope){
    if(addedIsotopes.includes(isotope)){//don't double add
        return null
    }
    else if(removedIsotopes.includes(isotope)){//avoid both adding and removing
        for(var i=0;i<removedIsotopes.length;i++){
            if(removedIsotopes[i] === isotope){
                removedIsotopes.splice(i,1);
                break;
            }
        }
    }
    else{
        addedIsotopes.push(isotope);
    }
    //add it to the list
    newElement = "<p class='iso-label'>" + isotope  + "</p>"
    removeButton = '<img class="rmv-btn" src="/icons/file-x.svg" onclick="removeIsotope('+"'"+isotope+"'"+')">'
    var ufl = document.getElementById("selectedIsotopes");
    ufl.innerHTML += DOMPurify.sanitize("<li class='list-group-item w-50' id='"+isotope+"'>" + newElement + removeButton +"</li>");
    document.getElementById("search-input").value = "";
    applyFilter();
}

/**
 * Remove an isotope from the analysis
 * @param {String} name 
 */
function removeIsotope(name){
    document.getElementById(name).remove();
    if(removedIsotopes.includes(name)){
        return null;
    }
    else if(addedIsotopes.includes(name)){
        for(var i=0;i<addedIsotopes.length;i++){
            if(addedIsotopes[i] === name){
                addedIsotopes.splice(i,1);
                break;
            }
        }
    }
    else{
        removedIsotopes.push(name);
    }
}

/**
 * Update which isotopes from the standard file are shown based on the search term
 */
function applyFilter(){
    var input = document.getElementById("search-input");
    var filter = input.value.toUpperCase();
    ul = document.getElementById("allIsotopes");
    li = ul.getElementsByTagName('li');
    for (i = 0; i < li.length; i++) {
        a = li[i].getElementsByClassName("iso-label")[0];
        txtValue = a.textContent || a.innerText;
        if (filter.length == 0 || !(txtValue.toUpperCase().startsWith(filter))) {//length = 0 is where no search has been entered, so no results should appear
        li[i].style.display = "none";//hide
        } else {
        li[i].style.display = "";//show
        }
    }
}

/**
 * Update the times in the NAA window when the file is changed
 */
function updateShownTimes(){
    var times = NAATimes[filesList.indexOf(document.getElementById("timeFileSelect").value)];
    if(times.length >= 1){
        document.getElementById("irrTimeInput").value = times[0].toString();
        document.getElementById("waitTimeInput").value = times[1].toString();
    }
}
