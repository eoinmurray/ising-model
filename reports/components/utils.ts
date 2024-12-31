

export function reshapeArray(array: any[], rows: number, columns: number) {
  if (rows * columns !== array.length) {
      throw new Error("The total elements do not match the new shape.");
  }
  
  let reshapedArray: any[][] = [];
  for (let i = 0; i < rows; i++) {
      reshapedArray.push(array.slice(i * columns, i * columns + columns));
  }
  return reshapedArray;
}

export function arrayToTidyData(array) {
  if (!Array.isArray(array) || array.length !== 10 || array.some(row => row.length !== 10)) {
      throw new Error("Input must be a 10x10 array.");
  }

  let tidyData: { x: number; y: number; value: any }[] = [];
  for (let rowIndex = 0; rowIndex < 10; rowIndex++) {
      for (let colIndex = 0; colIndex < 10; colIndex++) {
          tidyData.push({
              x: rowIndex,
              y: colIndex,
              value: array[rowIndex][colIndex]
          });
      }
  }
  return tidyData;
}
