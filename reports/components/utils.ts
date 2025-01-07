
type TidyRow = { row: number; column: number; value: number };

export function toTidyData(array: number[][]): TidyRow[] {
  const tidyData: TidyRow[] = [];
  
  array.forEach((row, rowIndex) => {
    row.forEach((value, colIndex) => {
      tidyData.push({ row: rowIndex, column: colIndex, value });
    });
  });

  return tidyData;
}

export function reshapeArray<T>(array: T[], rows: number, cols: number): T[][] {
  if (array.length !== rows * cols) {
    throw new Error("Array size does not match the specified dimensions.");
  }

  const reshaped: T[][] = [];
  for (let i = 0; i < rows; i++) {
    reshaped.push(array.slice(i * cols, i * cols + cols));
  }

  return reshaped;
}
