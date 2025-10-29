// const capsule1 = document.querySelector(".capsule-1");
// const capsule2 = document.querySelector(".capsule-2");
// const eigenValuesSection = document.querySelector(".eigenValues");
// const linearEquationsSection = document.querySelector(".linearEquations");
// const titleHeading = document.querySelector(".title");
// const submit1 = document.querySelector("#submit-1");
// const matrix1 = document.querySelectorAll('.eigenValues input[type="number"]');
// const displayMatrix = document.querySelectorAll(".displayMatrix .A");
// const outputDiv = document.querySelector(".output");
// const eigenValueDisplay = document.querySelectorAll(".eigenValue");


// let isActive = [true, false];
// let matrixA = [
//     [0, 0, 0],
//     [0, 0, 0],
//     [0, 0, 0]
// ];
// let eigenValues = [];
// let uniqueEigenValues = [];
// let repeatedEigenValues = [];
// let eigenValuesMap = new Map();
// const EPSILON = 1e-10; // Tolerance for floating-point comparison
// const eigenVectorsMap = new Map();


// capsule1.addEventListener("click", () => {
//     capsule1.classList.add("insetShadow");
//     capsule2.classList.remove("insetShadow");
//     titleHeading.textContent = "Eigen Values and Eigen Vector Calculator";
//     eigenValuesSection.classList.remove("hidden");
//     linearEquationsSection.classList.add("hidden");
// });

// capsule2.addEventListener("click", () => {
//     capsule2.classList.add("insetShadow");
//     capsule1.classList.remove("insetShadow");
//     titleHeading.textContent = "Linear Equations Solver";
//     linearEquationsSection.classList.remove("hidden");
//     eigenValuesSection.classList.add("hidden");
// });

// submit1.addEventListener(
//     "click",
//     calculateEigenValuesAndEigenVectors
// );














// function calculateEigenValuesAndEigenVectors() {
//     showOutput();
//     characteristicEquation();
//     countEigenValues();
//     findEigenVectors();
//     displayEigenVectors();
// }

// function displayEigenVectors() {

// }


// function solveLinearSystem(A, b) {
//     const n = A.length;

//     // 1. Create the augmented matrix [A | b]
//     const aug = new Array(n);
//     for (let i = 0; i < n; i++) {
//         aug[i] = [...A[i], b[i]];
//     }

//     // 2. Perform Gaussian Elimination (forward elimination)
//     for (let i = 0; i < n; i++) {
//         // Find pivot (largest absolute value in the column)
//         let maxRow = i;
//         for (let k = i + 1; k < n; k++) {
//             if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) {
//                 maxRow = k;
//             }
//         }

//         // Swap current row with pivot row
//         [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];

//         // Handle singular or near-singular matrix (should not happen for a valid chain)
//         if (Math.abs(aug[i][i]) < EPSILON) {
//             console.error("Matrix is singular, cannot find a unique solution.");
//             return [NaN, NaN, NaN]; // Return NaN if no solution
//         }

//         // Make pivot element 1
//         let pivot = aug[i][i];
//         for (let j = i; j <= n; j++) {
//             aug[i][j] /= pivot;
//         }

//         // Eliminate other rows
//         for (let k = 0; k < n; k++) {
//             if (k !== i) {
//                 let factor = aug[k][i];
//                 for (let j = i; j <= n; j++) {
//                     aug[k][j] -= factor * aug[i][j];
//                 }
//             }
//         }
//     }

//     // 3. Extract solution
//     // After Gauss-Jordan, the solution is the last column
//     const solution = new Array(n);
//     for (let i = 0; i < n; i++) {
//         solution[i] = aug[i][n];
//     }

//     return solution;
// }

// function findEigenVectors() {
//     eigenVectorsMap.clear();

//     // Iterate over unique eigenvalues and their counts (Algebraic Multiplicity)
//     for (const [lambda, am] of eigenValuesMap.entries()) {

//         // 1. Construct B = (A - λI)
//         const B = [
//             [matrixA[0][0] - lambda, matrixA[0][1], matrixA[0][2]],
//             [matrixA[1][0], matrixA[1][1] - lambda, matrixA[1][2]],
//             [matrixA[2][0], matrixA[2][1], matrixA[2][2]]
//         ];

//         // 2. Find the real eigenvectors (the null space of B)
//         // This gives us the Geometric Multiplicity (GM)
//         const basis = findEigenspaceBasis(matrixA, lambda); // Your existing function
//         const gm = basis.length;

//         // 3. Check if the matrix is defective
//         if (gm === am) {
//             // NOT defective. Just store the real eigenvectors.
//             eigenVectorsMap.set(lambda, basis);
//         } else {
//             // DEFECTIVE (GM < AM). We must build a Jordan Chain.
//             // This logic follows your Question2.java file.
//             console.log(`λ=${lambda} is defective (GM=${gm}, AM=${am}). Finding generalized vectors...`);

//             let chain = [...basis]; // Start the chain with the real eigenvector(s)

//             // We need to find (am - gm) more generalized vectors
//             // We'll build the chain off the first real eigenvector found
//             let v_prev = basis[0];

//             for (let i = gm; i < am; i++) {
//                 // Solve B * v_next = v_prev
//                 const v_next = solveLinearSystem(B, v_prev);

//                 // Add the new generalized vector to our chain
//                 chain.push(v_next);

//                 // This new vector becomes the right-hand side for the *next* iteration
//                 v_prev = v_next;
//             }

//             // Store the complete chain (real + generalized)
//             eigenVectorsMap.set(lambda, chain);
//         }
//     }

//     console.log("Eigenvectors Map:", eigenVectorsMap);
//     // Add your display logic here
// }

// // function findEigenVectors() {
// //     // Clear any previous results
// //     eigenVectorsMap.clear();

// //     // Iterate over the unique eigenvalues
// //     for (const lambda of eigenValuesMap.keys()) {
// //         // Find the basis (the set of eigenvectors) for this eigenvalue's eigenspace
// //         const basisVectors = findEigenspaceBasis(matrixA, lambda);

// //         // Store the array of eigenvectors
// //         eigenVectorsMap.set(lambda, basisVectors);
// //     }

// //     console.log("Eigenvectors Map:", eigenVectorsMap);
// //     // You can now add your display logic here
// //     // displayEigenVectors();
// // }

// function findEigenspaceBasis(A, lambda) {
//     // 1. Construct the matrix B = (A - λI)
//     const B = [
//         [A[0][0] - lambda, A[0][1], A[0][2]],
//         [A[1][0], A[1][1] - lambda, A[1][2]],
//         [A[2][0], A[2][1], A[2][2] - lambda]
//     ];

//     const r1 = B[0];
//     const r2 = B[1];
//     const r3 = B[2];

//     // 2. Try to find eigenvectors using the cross product of rows.
//     // A vector orthogonal to two rows is a solution.
//     const v1 = crossProduct(r1, r2);
//     const v2 = crossProduct(r1, r3);
//     const v3 = crossProduct(r2, r3);

//     const mag1 = magnitude(v1);
//     const mag2 = magnitude(v2);
//     const mag3 = magnitude(v3);

//     // 3. CASE 1: Geometric Multiplicity = 1 (Rank 2)
//     // The cross products are non-zero and parallel.
//     // We find the one with the largest magnitude for best precision.
//     if (mag1 > EPSILON) {
//         return [normalize(v1)]; // Found one eigenvector
//     }
//     if (mag2 > EPSILON) {
//         return [normalize(v2)]; // Found one eigenvector
//     }
//     if (mag3 > EPSILON) {
//         return [normalize(v3)]; // Found one eigenvector
//     }

//     // 4. CASE 2: Geometric Multiplicity > 1 (Rank 1 or 0)
//     // All cross products are zero, meaning all rows are parallel.

//     // Find the row with the largest magnitude (the most "significant" equation)
//     let r_main = r1;
//     if (magnitude(r2) > magnitude(r_main)) r_main = r2;
//     if (magnitude(r3) > magnitude(r_main)) r_main = r3;

//     // 4a. CASE 3: Geometric Multiplicity = 3 (Rank 0)
//     // The "most significant" row is also zero. (B is the zero matrix)
//     // This means A was A = λI. Any vector is an eigenvector.
//     if (magnitude(r_main) < EPSILON) {
//         return [
//             [1, 0, 0],
//             [0, 1, 0],
//             [0, 0, 1]
//         ];
//     }

//     // 4b. CASE 2 (cont'd): Geometric Multiplicity = 2 (Rank 1)
//     // We have one equation: x*r_main[0] + y*r_main[1] + z*r_main[2] = 0
//     // We need to find two independent vectors that satisfy this.
//     const [x, y, z] = r_main;
//     const absX = Math.abs(x);
//     const absY = Math.abs(y);
//     const absZ = Math.abs(z);

//     let basisVec1, basisVec2;

//     // Find the component with the largest absolute value to use as the pivot
//     if (absX >= absY && absX >= absZ) {
//         // Pivot on x: solve for x in terms of y and z
//         // Set y=1, z=0  => v1 = (-y/x, 1, 0)
//         basisVec1 = normalize([-y / x, 1, 0]);
//         // Set y=0, z=1  => v2 = (-z/x, 0, 1)
//         basisVec2 = normalize([-z / x, 0, 1]);
//     } else if (absY >= absX && absY >= absZ) {
//         // Pivot on y
//         // Set x=1, z=0  => v1 = (1, -x/y, 0)
//         basisVec1 = normalize([1, -x / y, 0]);
//         // Set x=0, z=1  => v2 = (0, -z/y, 1)
//         basisVec2 = normalize([0, -z / y, 1]);
//     } else {
//         // Pivot on z
//         // Set x=1, y=0  => v1 = (1, 0, -x/z)
//         basisVec1 = normalize([1, 0, -x / z]);
//         // Set x=0, y=1  => v2 = (0, 1, -y/z)
//         basisVec2 = normalize([0, 1, -y / z]);
//     }

//     return [basisVec1, basisVec2];
// }

// /**
//  * Calculates the cross product of two 3D vectors.
//  * v1 x v2
//  */
// function crossProduct(v1, v2) {
//     return [
//         v1[1] * v2[2] - v1[2] * v2[1], // x
//         v1[2] * v2[0] - v1[0] * v2[2], // y
//         v1[0] * v2[1] - v1[1] * v2[0]  // z
//     ];
// }

// /**
//  * Calculates the magnitude (L2 norm) of a vector.
//  */
// function magnitude(v) {
//     return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
// }

// /**
//  * Normalizes a vector (makes its length 1).
//  */
// function normalize(v) {
//     const mag = magnitude(v);
//     if (mag < EPSILON) {
//         return [0, 0, 0]; // Avoid division by zero
//     }
//     return [v[0] / mag, v[1] / mag, v[2] / mag];
// }
// function countEigenValues() {
//     // --- FIX 1: Use the global map, don't create a local 'countedMap' ---
//     eigenValuesMap.clear();
//     const sortedEigenValues = eigenValues.sort((a, b) => a - b);

//     uniqueEigenValues = [];
//     repeatedEigenValues = [];

//     // 1. Count Eigenvalues (Grouping numerically close values)
//     for (let ev of sortedEigenValues) {
//         let isNewValue = true;

//         // Iterate over existing keys to check for proximity
//         // --- FIX 2: Read keys from the global 'eigenValuesMap' ---
//         for (let existingKey of eigenValuesMap.keys()) {
//             if (Math.abs(ev - existingKey) < EPSILON) {
//                 // --- FIX 3: Set values on the global 'eigenValuesMap' ---
//                 eigenValuesMap.set(existingKey, eigenValuesMap.get(existingKey) + 1);
//                 isNewValue = false;
//                 break;
//             }
//         }

//         if (isNewValue) {
//             // --- FIX 4: Set values on the global 'eigenValuesMap' ---
//             eigenValuesMap.set(ev, 1);
//         }
//     }

//     // 2. Separate Unique and Repeated Eigenvalues
//     // --- FIX 5: Iterate over the global 'eigenValuesMap' ---
//     for (const [eigenvalue, count] of eigenValuesMap.entries()) {

//         // --- (See Potential Problem 2 below) ---
//         const roundedEv = Math.round(eigenvalue * 1e10) / 1e10; // Better rounding

//         if (count > 1) {
//             for (let i = 0; i < count; i++) {
//                 repeatedEigenValues.push(roundedEv);
//             }
//         } else {
//             uniqueEigenValues.push(roundedEv);
//         }
//     }

//     console.log("Eigen values map (Approx. Keys): ", eigenValuesMap); // Now logs the global map
//     console.log("Repeated eigen values: ", JSON.stringify(repeatedEigenValues));
//     console.log("Unique eigen values: ", uniqueEigenValues);
// }



// // function countEigenValues() {
// //     eigenValuesMap = new Map();
// //     for (let ev of eigenValues) {
// //         eigenValuesMap.set(ev, eigenValuesMap.get(ev) == undefined ? 1 : eigenValuesMap.get(ev) + 1);
// //     }

// //     for (let k of eigenValuesMap) {
// //         let count = eigenValuesMap.get(k);
// //         if (count > 1) {
// //             if (count == 2) repeatedEigenValues = [k, k];
// //             else repeatedEigenValues = [k, k, k]
// //         } else {
// //             uniqueEigenValues.push(k);
// //         }
// //     }
// //     console.log("Eigen values map : ", eigenValuesMap);
// //     console.log("repeated eigen values ", JSON.stringify(repeatedEigenValues));
// //     console.log("Unique eigen values : ", uniqueEigenValues);
// // }

// function characteristicEquation() {
//     let traceA = matrixA[0][0] + matrixA[1][1] + matrixA[2][2];
//     let detA = determinant(matrixA);
//     let sumOfMinors = minorSum(matrixA);
//     solveCubic(1, -traceA, sumOfMinors, -detA);
//     console.log([1, traceA, sumOfMinors, detA]);
// }

// // function solveCubic(a, b, c, d) {
// //     let p = (3 * a * c - Math.pow(b, 2)) / (3 * Math.pow(a, 2));
// //     let q = (2 * Math.pow(b, 3) - 9 * a * b * c + 27 * Math.pow(a, 2) * d) / (27 * Math.pow(a, 3));
// //     let delta = Math.pow(q / 2, 2) + Math.pow(p / 3, 3);
// //     console.log(delta);
// //     const offset = b / (3 * a);
// //     findRoots(p, q, delta, offset);
// //     displayEigenValues();

// // }

// function displayEigenValues() {
//     eigenValues.sort();
//     let i = 0;
//     for (let ev of eigenValueDisplay) {
//         ev.textContent = Math.round(eigenValues[i++]);
//     }
// }

// function solveCubic(a, b, c, d) {
//     let p = (3 * a * c - Math.pow(b, 2)) / (3 * Math.pow(a, 2));
//     let q = (2 * Math.pow(b, 3) - 9 * a * b * c + 27 * Math.pow(a, 2) * d) / (27 * Math.pow(a, 3));
//     const offset = b / (3 * a);

//     // Define a small tolerance for zero
//     const EPSILON = 1e-10;

//     // --- NEW FIX: Handle the special case where q is near zero ---
//     if (Math.abs(q) < EPSILON) {
//         console.log("case (q=0 special)");
//         // Equation is y^3 + py = 0
//         // Roots are 0, sqrt(-p), -sqrt(-p)
//         const y1 = 0;
//         const y2 = Math.sqrt(-p); // p must be negative if original roots are real
//         const y3 = -Math.sqrt(-p);

//         setEigenValues(y1, y2, y3, offset);
//         displayEigenValues(); // Make sure to call your display function
//         return; // Exit function early
//     }
//     // --- End of new fix ---

//     let delta = Math.pow(q / 2, 2) + Math.pow(p / 3, 3);
//     console.log(delta);

//     // We still need the delta check for the *other* case (repeated roots)
//     if (Math.abs(delta) < EPSILON) {
//         delta = 0; // Force it to zero
//     }

//     findRoots(p, q, delta, offset);
//     displayEigenValues();
// }

// // Your findRoots function (with the delta === 0 check) is now fine
// function findRoots(p, q, delta, offset) {
//     if (delta > 0) {
//         console.log("case 1");
//         const u = Math.cbrt(-q / 2 + Math.sqrt(delta));
//         const v = Math.cbrt(-q / 2 - Math.sqrt(delta));
//         const y1 = u + v;

//         // This case has 1 real and 2 complex roots
//         setEigenValues(y1, NaN, NaN, offset);
//     } else if (delta === 0) { // This will now catch the tiny delta
//         console.log("case 2");
//         const y_root = -Math.cbrt(-q / 2);
//         // Roots are y_root (repeated) and -2*y_root (single)
//         setEigenValues(y_root, y_root, -2 * y_root, offset);
//     } else {
//         // Case 3 (delta < 0)
//         // This case will no longer be hit by (q=0), so it's more stable
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

// // function findRoots(p, q, delta, offset) {
// //     // Define a small tolerance for zero
// //     const EPSILON = 1e-10;

// //     // --- FIX 1: Check if delta is CLOSE to zero ---
// //     if (Math.abs(delta) < EPSILON) {
// //         console.log("case 2 (Delta is effectively zero)");

// //         const y1 = Math.cbrt(-q / 2) * 2;
// //         const y2 = -Math.cbrt(-q / 2);

// //         // Pass the roots: one root (y1) and the repeated root (y2)
// //         // Your setEigenValues expects three, so we send the repeated root twice.
// //         setEigenValues(y2, y2, y1, offset);

// //         // --- FIX 2: Check if delta is clearly POSITIVE ---
// //     } else if (delta > 0) {
// //         console.log("case 1");
// //         const u = Math.cbrt(-q / 2 + Math.sqrt(delta));
// //         const v = Math.cbrt(-q / 2 - Math.sqrt(delta));
// //         const y1 = u + v;

// //         // --- FIX 3: Your logic here was also incorrect ---
// //         // u and v are NOT roots, they are parts of the calculation.
// //         // In Case 1, there is only ONE real root (y1). 
// //         // The other two roots are complex. Your code was incorrectly 
// //         // passing u and v as if they were roots.

// //         // We'll pass y1 as the real root. The other two 
// //         // will be NaN (Not a Number) since they are complex.
// //         setEigenValues(y1, NaN, NaN, offset);

// //     } else {
// //         // Case 3 (delta < 0)
// //         console.log("case 3");

// //         const TWO_PI = 2 * Math.PI;
// //         const r = Math.sqrt(-Math.pow(p / 3, 3));
// //         const ratio = (-q / 2) / r;
// //         const clampedRatio = Math.min(1, Math.max(-1, ratio));
// //         const phi = Math.acos(clampedRatio);
// //         const factor = 2 * Math.cbrt(r);

// //         const y1 = factor * Math.cos(phi / 3);
// //         const y2 = factor * Math.cos((phi + TWO_PI) / 3);
// //         const y3 = factor * Math.cos((phi + 4 * Math.PI) / 3);

// //         setEigenValues(y1, y2, y3, offset);
// //     }
// // }


// // function findRoots(p, q, delta, offset) {
// //     if (delta > 0) {
// //         console.log("case 1");
// //         const u = Math.cbrt(-q / 2 + Math.sqrt(delta));
// //         const v = Math.cbrt(-q / 2 - Math.sqrt(delta));
// //         const y1 = u + v;

// //         setEigenValues(u, v, y1, offset);
// //     } else if (delta === 0) {
// //         console.log("case 2");

// //         const y1 = Math.cbrt(-q / 2) * 2;
// //         const y2 = -Math.cbrt(-q / 2);

// //         setEigenValues(y1, y1, y2, offset);
// //     } else {
// //         console.log("case 3");

// //         const TWO_PI = 2 * Math.PI;

// //         const r = Math.sqrt(-Math.pow(p / 3, 3));

// //         const ratio = (-q / 2) / r;
// //         const clampedRatio = Math.min(1, Math.max(-1, ratio));
// //         const phi = Math.acos(clampedRatio);

// //         const factor = 2 * Math.cbrt(r);

// //         const y1 = factor * Math.cos(phi / 3);
// //         const y2 = factor * Math.cos((phi + TWO_PI) / 3);
// //         const y3 = factor * Math.cos((phi + 4 * Math.PI) / 3);

// //         setEigenValues(y1, y2, y3, offset);
// //     }
// // }

// function setEigenValues(y1, y2, y3, offset) {
//     eigenValues = [y1 - offset, y2 - offset, y3 - offset];
//     console.log("Eigen values are : ", eigenValues);
// }

// function minorSum(x) {
//     // $c = (a_{22}a_{33} - a_{23}a_{32}) + (a_{11}a_{33} - a_{13}a_{31}) + (a_{11}a_{22} - a_{12}a_{21})$
//     return (x[1][1] * x[2][2] - x[1][2] * x[2][1]) + (x[0][0] * x[2][2] - x[0][2] * x[2][0]) + (x[0][0] * x[1][1] - x[0][1] * x[1][0]);
// }

// function determinant(x) {
//     let one = x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1]);
//     let two = -x[0][1] * (x[1][0] * x[2][2] - x[1][2] * x[2][0]);
//     let three = x[0][2] * (x[1][0] * x[2][1] - x[2][0] * x[1][1]);
//     return one + two + three;
// }

// function showOutput() {
//     let i = 0;
//     outputDiv.hidden = false;
//     for (let inp of matrix1) {
//         matrixA[Math.floor(i / 3)][i % 3] = parseInt(inp.value | 0);
//         displayMatrix[i].textContent = matrixA[Math.floor(i / 3)][i % 3];
//         i++;
//     }
//     console.log(matrixA);
// }

