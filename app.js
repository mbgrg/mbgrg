var Module = {
    preRun: [],
    postRun: [],
    print: (function () {
        let element = document.getElementById("output");
        if (element) element.value = ""; // clear browser cache
        return function (text) {
            if (arguments.length > 1)
                text = Array.prototype.slice.call(arguments).join(" ");
            if (element) {
                element.value += text + "\n";
                element.scrollTop = element.scrollHeight; // focus on bottom
            }
        };
    })(),
    canvas: (() => {
        return document.getElementById("canvas");
    })(),
    setStatus: (text) => {
        if (!Module.setStatus.last)
            Module.setStatus.last = { time: Date.now(), text: "" };
        if (text === Module.setStatus.last.text) return;
        let m = text.match(/([^(]+)\((\d+(\.\d+)?)\/(\d+)\)/);
        let now = Date.now();
        if (m && now - Module.setStatus.last.time < 30) return; // if this is a progress update, skip it if too soon
        Module.setStatus.last.time = now;
        Module.setStatus.last.text = text;
        if (m) {
            text = m[1];
            progressElement.value = parseInt(m[2]) * 100;
            progressElement.max = parseInt(m[4]) * 100;
            progressElement.hidden = false;
            spinnerElement.hidden = false;
        } else {
            progressElement.value = null;
            progressElement.max = null;
            progressElement.hidden = true;
            if (!text) spinnerElement.hidden = true;
        }
        statusElement.innerHTML = text;
    },
    totalDependencies: 0,
    monitorRunDependencies: (left) => {
        this.totalDependencies = Math.max(this.totalDependencies, left);
        Module.setStatus(
            left
                ? "Preparing... (" +
                (this.totalDependencies - left) +
                "/" +
                this.totalDependencies +
                ")"
                : "All downloads complete.",
        );
    },
};
Module.setStatus("Downloading...");
window.onerror = (message, source, lineno, colno, error) => {
    if (
        message.includes('Uncaught TypeError: Module.ccall is not a function') || message.includes('Uncaught RuntimeError: Aborted(Assertion failed: native function `stackSave` called before runtime initialization)')
    ) {
        // Espera 5 segundos antes de recargar la página
        setTimeout(() => {
            console.log("hola");
            window.location.reload();
        }, 500);
    } else {
        // Si no es el error específico, solo muestra el mensaje de error en la consola
        console.error(message);
    }

    Module.setStatus("Exception thrown, see JavaScript console: " + message);
    spinnerElement.style.display = "none";
    Module.setStatus = (text) => {
        if (text) console.error("[post-exception status] " + text);
    };
};


function computeShape() {
    svgOutput.value = Module.ccall(
        "_Z19ComputeSVGFromShapePc", // name of C function
        "string", // return type
        //['number','string'], // argument types
        ["string"], // argument types
        [shapeInput.value], // arguments
    );
    drawFigure();
}
// para cargar el shape.txt en el textarea al principio
window.addEventListener("load", function () {
    // Crear una nueva solicitud XMLHttpRequest
    let solicitud = new XMLHttpRequest();

    // Configurar la solicitud
    solicitud.open("GET", "./shape.txt", true);
    solicitud.onreadystatechange = function () {
        if (solicitud.readyState === 4 && solicitud.status === 200) {
            // Obtener el contenido del archivo
            let contenido = solicitud.responseText;

            // Establecer el contenido en el textarea
            let textarea = document.getElementById("shapeInput");
            textarea.value = contenido;
        }
    };

    // Enviar la solicitud
    solicitud.send();
});


function erasePoint(i) {
    shapeInput.value = Module.ccall("_Z10erasepointi", "string", ["number"], [i]);
    computeShape();
}

function preview() {
    return Module.ccall("_Z11downloadsvgv", "string", ["number"], []);
}

function downloadsvg() {
    var aux = Module.ccall("_Z11downloadsvgv", "string", ["number"], []);

    var enlaceDescarga = document.createElement("a");
    enlaceDescarga.setAttribute(
        "href",
        "data:image/svg+xml;charset=utf-8," + encodeURIComponent(aux),
    );
    enlaceDescarga.setAttribute("download", "shape.svg");
    enlaceDescarga.style.display = "none";

    document.body.appendChild(enlaceDescarga);
    enlaceDescarga.click();
    document.body.removeChild(enlaceDescarga);
}

function randomGenerate() {
    handleNoActionClick();
    // Antes de generar un nuevo SVG, verifica si ya hay uno previo y libera la memoria
    if (typeof svgPtr !== "undefined") {
        Module["_free"](svgPtr); // Liberar la memoria del SVG anterior
    }

    shapeInput.value = Module.ccall("_Z16randomgenerationv", "string");
    computeShape();

    // Guarda el nuevo puntero a la memoria del SVG generado
    svgPtr = Module.ccall(
        "_Z19ComputeSVGFromShapePc",
        "number",
        ["string"],
        [shapeInput.value],
    );
    fixit();
    addShape();
}

function similarGenerate() {
    handleNoActionClick();
    // Antes de generar un nuevo SVG, verifica si ya hay uno previo y libera la memoria
    if (typeof svgPtr !== "undefined") {
        Module["_free"](svgPtr); // Liberar la memoria del SVG anterior
    }

    shapeInput.value = Module.ccall("_Z17similargenerationv", "string");
    computeShape();

    // Guarda el nuevo puntero a la memoria del SVG generado
    svgPtr = Module.ccall(
        "_Z19ComputeSVGFromShapePc",
        "number",
        ["string"],
        [shapeInput.value],
    );
    fixit();
    addShape();
}

function createconection(circleIndex1, circleIndex2) {
    shapeInput.value = Module.ccall(
        "_Z12connectnodesii",
        "string",
        ["number", "number"],
        [circleIndex1, circleIndex2],
    );
    computeShape();
}

function eraseConnection(circleIndex1, circleIndex2) {
    shapeInput.value = Module.ccall(
        "_Z15disconnectnodesdd",
        "string",
        ["number", "number"],
        [circleIndex1, circleIndex2],
    );
    computeShape();
}

function middleCircle(x, y) {
    shapeInput.value = Module.ccall(
        "_Z17insertpointmiddledd",
        "string",
        ["number", "number"],
        [x, y],
    );
    computeShape();
}

function addpoint(x, r) {
    shapeInput.value = Module.ccall(
        "_Z11insertpointii",
        "string",
        ["number", "number"],
        [x, r],
    );
    computeShape();
}

function transform(zoom_factor, zx, zy, dx, dy) {
    shapeInput.value = Module.ccall(
        "_Z9transformfffff",
        "string",
        ["float", "float", "float", "float", "float"],
        [zoom_factor, zx, zy, dx, dy],
    );
}

window.onload = function() {
    // Retrasa la simulación del clic por 100 milisegundos (ajusta según sea necesario)
    setTimeout(function() {
        // Obtén el elemento del botón por su ID
        var computeButton = document.getElementById("computeButton");

        // Simula un clic en el botón
        computeButton.click();
    }, 100);
};