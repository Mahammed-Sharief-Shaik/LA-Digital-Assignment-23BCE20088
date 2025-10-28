const capsule1 = document.querySelector(".capsule-1");
const capsule2 = document.querySelector(".capsule-2");
const eigenValuesSection = document.querySelector(".eigenValues");
const linearEquationsSection = document.querySelector(".linearEquations");
const titleHeading = document.querySelector(".title");
const submit1 = document.querySelector("#submit-1");
const matrix1 = document.querySelectorAll('.eigenValues input[type="number"]');
const displayMatrix = document.querySelectorAll(".displayMatrix .A");
const outputDiv = document.querySelector(".output");



let isActive = [true, false];
let matrixA = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
];


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
    let i = 0;
    outputDiv.hidden = false;
    for (let inp of matrix1) {
        matrixA[Math.floor(i / 3)][i % 3] = parseInt(inp.value | 0);
        displayMatrix[i].textContent = matrixA[Math.floor(i / 3)][i % 3];
        i++;
    }

    console.log(matrixA);

}