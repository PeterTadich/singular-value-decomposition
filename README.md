# singular-value-decomposition
Matrix decomposition.

## Dependencies

There is 1 dependency 'matrix-computations'.

```bash
https://github.com/PeterTadich/matrix-computations
```

## Installation

### Node.js

```bash
npm install https://github.com/PeterTadich/singular-value-decomposition#main
```

### Google Chrome Web browser

No installation required for the Google Chrome Web browser.

## How to use

### Node.js/Google Chrome Web browser

```js
import * as svdcmp from 'singular-value-decomposition';
```

## Examples

### Node.js (server side)

Copy the following code to index.mjs

```js
import * as svdcmp from 'singular-value-decomposition';
import * as hlao from 'matrix-computations';

var a = [[0,0,0,0],[0,2,3,5],[0,-4,2,3]]; //where row index = 0 is 0 or 'undefined' and column index 0 is 0 or 'undefined'
var m = 2; //number of rows
var n = 3; //number of columns
var w = hlao.zeros_vector((n+1),'row'); //row vector where index = 0 is undefined
var v = hlao.zeros_matrix((n+1),(n+1)); //matrix
var uwv = svdcmp.svdcmp(a, m, n, w, v);
console.log(uwv);
```

Then run:

```bash
npm init -y
npm install https://github.com/PeterTadich/singular-value-decomposition#main https://github.com/PeterTadich/matrix-computations
node index.mjs
```

If the above does not work, modify the package.json file as follows:
Helpful ref: [https://stackoverflow.com/questions/45854169/how-can-i-use-an-es6-import-in-node-js](https://stackoverflow.com/questions/45854169/how-can-i-use-an-es6-import-in-node-js)

```js
"scripts": {
    "test": "echo \"Error: no test specified\" && exit 1",
    "start": "node --experimental-modules index.mjs"
  },
"type": "module",
```

Now try:

```bash
npm start
```

Returns:

```js
[
  //u
  [
    [  0.0000,  0.0000,  0.0000,  0.0000 ],
    [  0.0000,  0.8146,  0.5800,  0.0000 ],
    [  0.0000,  0.5800, -0.8146,  0.0000 ]
  ],
  //w
  [  0.0000,  6.8744,  4.4433,  -0.0000 ],
  //v
  [
    [  0.0000,  0.0000,  0.0000,  0.0000 ],
    [  0.0000, -0.1005,  0.9944,  0.0327 ],
    [  0.0000,  0.5242,  0.0250,  0.8512 ],
    [  0.0000,  0.8456,  0.1027, -0.5238 ]
  ]
]
```

Notice: the first row/column need to be removed.

### Google Chrome Web browser (client side)

Copy the following code to index.mjs

```js
import * as svdcmp from './svdcmp.mjs';
import * as hlao from './hlao.mjs';

var a = [[0,0,0,0],[0,2,3,5],[0,-4,2,3]]; //where row index = 0 is 0 or 'undefined' and column index 0 is 0 or 'undefined'
var m = 2; //number of rows
var n = 3; //number of columns
var w = hlao.zeros_vector((n+1),'row'); //row vector where index = 0 is undefined
var v = hlao.zeros_matrix((n+1),(n+1)); //matrix
var uwv = svdcmp.svdcmp(a, m, n, w, v);
console.log(uwv);
```

Copy the following code to index.html

```html
<!DOCTYPE html>
<html lang="en">
    <head>
        <script type="module" src="index.mjs"></script>
    </head>
</html>
```

Download the modules:

```bash
curl -O https://raw.githubusercontent.com/PeterTadich/singular-value-decomposition/main/svdcmp.mjs
curl -O https://raw.githubusercontent.com/PeterTadich/matrix-computations/main/hlao.mjs
```

Start the server, open index.html with Google Chrome and inspect with the console.

Returns:

```js
[
  //u
  [
    [  0.0000,  0.0000,  0.0000,  0.0000 ],
    [  0.0000,  0.8146,  0.5800,  0.0000 ],
    [  0.0000,  0.5800, -0.8146,  0.0000 ]
  ],
  //w
  [  0.0000,  6.8744,  4.4433,  -0.0000 ],
  //v
  [
    [  0.0000,  0.0000,  0.0000,  0.0000 ],
    [  0.0000, -0.1005,  0.9944,  0.0327 ],
    [  0.0000,  0.5242,  0.0250,  0.8512 ],
    [  0.0000,  0.8456,  0.1027, -0.5238 ]
  ]
]
```

Server '*.conf' files may require modification.

mime.types file include:

```conf
text/javascript mjs js
```

httpd.conf file include:

```conf
AddType text/javascript .mjs
```

Alternatively change the file extension from *.mjs to *.js

## License

[MIT](LICENSE)