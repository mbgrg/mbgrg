//Event Handlers
function HandleExportShapeButtonClick() {
    const shapeContent = shapeInput.value;

    const blob = new Blob([shapeContent], {
        type: "text/plain",
    });
    const url = URL.createObjectURL(blob);

    const a = document.createElement("a");
    a.href = url;
    a.download = "shape.txt";
    a.style.display = "none";
    document.body.appendChild(a);

    a.click();

    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

function handleNoActionClick() {
    removeAllCanvasEvents();
    addAllCanvasEvents();
    setButtonStyle();
}
function handleBPreviewMouseDown() {
    svgOutput.value = preview();
    drawFigure();
}

function handleBPreviewMouseUp() {
    computeShape();
}

function handleMiddleCircleClick() {
    handleNoActionClick();
    buttonColor(middlecircle, "#ff9800");
    waitingMessage.style.display = "block";
    removeAllCanvasEvents();
    canvas.addEventListener("click", handleMiddleCircleConnectionClick);
}

function handleMiddleCircleConnectionClick() {
    var coords = relMouseCoords(event);

    middleCircle(coords.x , coords.y);
    waitingMessage.style.display = "none";
    buttonColor(middlecircle, "#4CAF50");
    canvas.removeEventListener("click", handleMiddleCircleConnectionClick);
    addShape();

}
function handleEraseconClick() {
    handleNoActionClick();
    buttonColor(erasecon, "#ff9800");
    waitingMessage.style.display = "block";
    removeAllCanvasEvents();
    canvas.addEventListener("click", handleEraseConnectionClick);
}

function handleEraseConnectionClick(event) {
    var coords = relMouseCoords(event);
    eraseConnection(coords.x, coords.y );
    waitingMessage.style.display = "none";
    buttonColor(erasecon, "#4CAF50");
    canvas.removeEventListener("click", handleEraseConnectionClick);
    addShape();
}
function handleCreateconClick() {
    handleNoActionClick();
    buttonColor(createcon, "#ff9800");
    if (circulitos) {
        waitingMessage.style.display = "block";
        removeAllCanvasEvents();
        circulitos.on("click", handleCreateConnectionClick);
    }
}

function handleCreateConnectionClick() {
    const circleIndex = Array.from(circulitos.nodes()).indexOf(this);

    if (typeof createconection.firstCircleIndex === "undefined") {
        createconection.firstCircleIndex = circleIndex;
    } else {
        createconection(createconection.firstCircleIndex, circleIndex);
        delete createconection.firstCircleIndex;
        waitingMessage.style.display = "none";
        buttonColor(createcon, "#4CAF50");
    }
    addShape();
}
function handleAddButtonClick() {
    handleNoActionClick();
    buttonColor(addButton, "#ff9800");
    waitingMessage.style.display = "block";
    removeAllCanvasEvents();
    canvas.addEventListener("click", addcircle, { once: true });
}

function handleEraseButtonClick() {
    handleNoActionClick();
    buttonColor(eraseButton, "#ff9800");
    if (circulitos) {
        waitingMessage.style.display = "block";
        removeAllCanvasEvents();
        circulitos.on("click", handleCircleClickForErase);
    }
}
function handleCircleClickForErase() {
    const circleIndex = Array.from(circulitos.nodes()).indexOf(this);
    erasePoint(circleIndex);
    waitingMessage.style.display = "none";
    buttonColor(eraseButton, "#4CAF50");
    addShape();
}
function handleComputeButtonClick() {
    addShape(shapeInput.value.split("\n"));
    computeShape();
    fixit();
}

function handleHelpButtonClick() {
    popuphelp.style.display = "block";
}
function handleCloseButtonClick() {
    popuphelp.style.display = "none";
}

function handleFileUpload(event) {
    const file = event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = function (event) {
        const fileContent = event.target.result;
        shapeInput.value = fileContent.trim();
        computeShape();
    };

    reader.readAsText(file, "UTF-8");
}
function handleWheelEvent(event) {
    if (isMouseOverSVG(event)) {
        event.preventDefault(); // Evita el desplazamiento predeterminado de la página
    }
}
function handleCanvasMouseDown(event) {
    if (!isRightClick(event) || popup.style.display === "block") return;
    var coords = relMouseCoords(event);
    initialClickX = coords.x;
    initialClickY = coords.y;
    isCanvasClicked = true;
}

function handleCanvasMouseUp() {
    isCanvasClicked = false;
    addShape();
}
let wheelTimer;

// Función para guardar la figura
function saveFigure() {
    addShape(); // Guardar la figura
    // Restablecer el temporizador
    clearTimeout(wheelTimer);
    wheelTimer = null;
}

// Controlador de eventos para la rueda del ratón en el canvas
function handleWheelCanvas(event) {
    if (selectedcircle == null) {
        var isScrollUp = event.deltaY < 0;
        var coords = relMouseCoords(event);
        var zoom = isScrollUp ? 1.03 : 0.97;
        transform(zoom, coords.x, coords.y, 0, 0);
        computeShape();

        // Reiniciar el temporizador cada vez que se usa la rueda del ratón
        clearTimeout(wheelTimer);
        wheelTimer = setTimeout(saveFigure, 300); // 300 milisegundos después de terminar el uso de la rueda
    }
}
function HandleMovecanvas(event) {
    if (!isRightClick(event) || popup.style.display === "block") return;
    if (isCanvasClicked) {
        if (activeCircle === false) {
            var coords = relMouseCoords(event);
            var currentMouseX = coords.x;
            var currentMouseY = coords.y;
            var movementX = currentMouseX - initialMouseX;
            var movementY = currentMouseY - initialMouseY;

            transform(1, null, null, movementX, movementY);
            computeShape();

            initialMouseX = currentMouseX;
            initialMouseY = currentMouseY;
        }
    }
}
function handleCanvasContextMenu(event) {
    event.preventDefault();
    openPopup(event.clientX, event.clientY);
}

function handleCanvasMouseDown2(event) {
    var coords = relMouseCoords(event);
    initialMouseX = coords.x;
    initialMouseY = coords.y;
}

function handleDocumentKeyDown(event) {
    if (event.ctrlKey && event.keyCode === 90) {
        updateShapeInput();
    } else if (event.ctrlKey && event.keyCode === 89) {
        redoShapeInput();
    }
}
