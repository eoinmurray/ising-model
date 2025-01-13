# Browse Simulation Runs

```ts
const dirs = await FileAttachment('list.json').json()

display(
  Inputs.table(
    dirs.directories, {
      format: {
        path: (x) => html`<a href="simulation/${x}">${x}</a>`
      }
    }
  )
)
```
