const capsule1 = document.querySelector(".capsule-1");
const capsule2 = document.querySelector(".capsule-2");
const eigenValuesSection = document.querySelector(".eigenValues");
const linearEquationsSection = document.querySelector(".linearEquations");
const titleHeading = document.querySelector(".title");
const submit1 = document.querySelector("#submit-1");
const matrix1 = document.querySelectorAll('.eigenValues input[type="number"]');
const displayMatrix = document.querySelectorAll(".displayMatrix .A");
const outputDiv = document.querySelector(".output");
const eigenValueDisplay = document.querySelectorAll(".eigenValue");


let isActive = [true, false];
let matrixA = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
];
let eigenValues = new Array(3);

capsule1.addEventListener("click", () => {
    capsule1.classList.add("insetShadow");
    capsule2.classList.remove("insetShadow");
    titleHeading.textContent = "Eigen Values and Eigen Vector Calculator";
    eigenValuesSection.classList.remove("hidden");
    linearEquationsSection.classList.add("hidden");
});

capsule2.addEventListener("click", () => {
    capsule2.classList.add("insetShadow");
    capsule1.classList.remove("insetShadow");
    titleHeading.textContent = "Linear Equations Solver";
    linearEquationsSection.classList.remove("hidden");
    eigenValuesSection.classList.add("hidden");
});

submit1.addEventListener(
    "click",
    calculateEigenValues
);














function calculateEigenValues() {
    showOutput();
    characteristicEquation();

}

function characteristicEquation() {
    let traceA = matrixA[0][0] + matrixA[1][1] + matrixA[2][2];
    let detA = determinant(matrixA);
    let sumOfMinors = minorSum(matrixA);
    solveCubic(1, -traceA, sumOfMinors, -detA);
    console.log([1, traceA, sumOfMinors, detA]);
}

function solveCubic(a, b, c, d) {
    let p = (3 * a * c - Math.pow(b, 2)) / (3 * Math.pow(a, 2));
    let q = (2 * Math.pow(b, 3) - 9 * a * b * c + 27 * Math.pow(a, 2) * d) / (27 * Math.pow(a, 3));
    let delta = Math.pow(q / 2, 2) + Math.pow(p / 3, 3);
    console.log(delta);
    const offset = b / (3 * a);
    findRoots(p, q, delta, offset);
    displayEigenValues();

}

function displayEigenValues() {
    let i = 0;
    for (let ev of eigenValueDisplay) {
        ev.textContent = eigenValues[i++];
    }
}

function findRoots(p, q, delta, offset) {
    // Define a small tolerance for zero
    const EPSILON = 1e-10;

    // --- FIX 1: Check if delta is CLOSE to zero ---
    if (Math.abs(delta) < EPSILON) {
        console.log("case 2 (Delta is effectively zero)");

        const y1 = Math.cbrt(-q / 2) * 2;
        const y2 = -Math.cbrt(-q / 2);

        // Pass the roots: one root (y1) and the repeated root (y2)
        // Your setEigenValues expects three, so we send the repeated root twice.
        setEigenValues(y2, y2, y1, offset);

        // --- FIX 2: Check if delta is clearly POSITIVE ---
    } else if (delta > 0) {
        console.log("case 1");
        const u = Math.cbrt(-q / 2 + Math.sqrt(delta));
        const v = Math.cbrt(-q / 2 - Math.sqrt(delta));
        const y1 = u + v;

        // --- FIX 3: Your logic here was also incorrect ---
        // u and v are NOT roots, they are parts of the calculation.
        // In Case 1, there is only ONE real root (y1). 
        // The other two roots are complex. Your code was incorrectly 
        // passing u and v as if they were roots.

        // We'll pass y1 as the real root. The other two 
        // will be NaN (Not a Number) since they are complex.
        setEigenValues(y1, NaN, NaN, offset);

    } else {
        // Case 3 (delta < 0)
        console.log("case 3");

        const TWO_PI = 2 * Math.PI;
        const r = Math.sqrt(-Math.pow(p / 3, 3));
        const ratio = (-q / 2) / r;
        const clampedRatio = Math.min(1, Math.max(-1, ratio));
        const phi = Math.acos(clampedRatio);
        const factor = 2 * Math.cbrt(r);

        const y1 = factor * Math.cos(phi / 3);
        const y2 = factor * Math.cos((phi + TWO_PI) / 3);
        const y3 = factor * Math.cos((phi + 4 * Math.PI) / 3);

        setEigenValues(y1, y2, y3, offset);
    }
}
// function findRoots(p, q, delta, offset) {
//     if (delta > 0) {
//         console.log("case 1");
//         const u = Math.cbrt(-q / 2 + Math.sqrt(delta));
//         const v = Math.cbrt(-q / 2 - Math.sqrt(delta));
//         const y1 = u + v;

//         setEigenValues(u, v, y1, offset);
//     } else if (delta === 0) {
//         console.log("case 2");

//         const y1 = Math.cbrt(-q / 2) * 2;
//         const y2 = -Math.cbrt(-q / 2);

//         setEigenValues(y1, y1, y2, offset);
//     } else {
//         console.log("case 3");

//         const TWO_PI = 2 * Math.PI;

//         const r = Math.sqrt(-Math.pow(p / 3, 3));

//         const ratio = (-q / 2) / r;
//         const clampedRatio = Math.min(1, Math.max(-1, ratio));
//         const phi = Math.acos(clampedRatio);

//         const factor = 2 * Math.cbrt(r);

//         const y1 = factor * Math.cos(phi / 3);
//         const y2 = factor * Math.cos((phi + TWO_PI) / 3);
//         const y3 = factor * Math.cos((phi + 4 * Math.PI) / 3);

//         setEigenValues(y1, y2, y3, offset);
//     }
// }

function setEigenValues(y1, y2, y3, offset) {
    eigenValues = [y1 - offset, y2 - offset, y3 - offset];
    console.log("Eigen values are : ", eigenValues);
}

function minorSum(x) {
    // $c = (a_{22}a_{33} - a_{23}a_{32}) + (a_{11}a_{33} - a_{13}a_{31}) + (a_{11}a_{22} - a_{12}a_{21})$
    return (x[1][1] * x[2][2] - x[1][2] * x[2][1]) + (x[0][0] * x[2][2] - x[0][2] * x[2][0]) + (x[0][0] * x[1][1] - x[0][1] * x[1][0]);
}

function determinant(x) {
    let one = x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1]);
    let two = -x[0][1] * (x[1][0] * x[2][2] - x[1][2] * x[2][0]);
    let three = x[0][2] * (x[1][0] * x[2][1] - x[2][0] * x[1][1]);

    return one + two + three;
}

function showOutput() {
    let i = 0;
    outputDiv.hidden = false;
    for (let inp of matrix1) {
        matrixA[Math.floor(i / 3)][i % 3] = parseInt(inp.value | 0);
        displayMatrix[i].textContent = matrixA[Math.floor(i / 3)][i % 3];
        i++;
    }

    console.log(matrixA);
}