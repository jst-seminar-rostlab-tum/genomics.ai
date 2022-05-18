const options = require("./emails.config.js");
const juice = require("juice");
const fs = require("fs/promises");
const path = require("path");

const tailwindcssfile = path.resolve(__dirname, "build/tailwind.css");
const tailwindcss = fs.readFile(tailwindcssfile);

(async function main() {
  console.log(options.input);
  const basePath = path.resolve(__dirname, options.input);
  const outputBase = path.resolve(__dirname, options.output);

  let dirs = await fs.readdir(basePath);
  console.log(dirs);
  for (const dir of dirs) {
    await processFolder(path.resolve(basePath, dir), path.resolve(outputBase, dir));
  }
})();

async function processFolder(basePath, outputBase) {
  let files = await listFiles(path.resolve(basePath));
  let htmls = files.filter((f) => f.endsWith(".html"));
  let others = files.filter((f) => !f.endsWith(".html"));
  for (const html of htmls) {
    await compileEmail(basePath, html, outputBase);
  }
  console.log("Processing "+basePath);
  for (const other of others) {
    console.log("Copy "+other);
    let parentdirs = path.dirname(other);
    await fs.mkdir(path.resolve(outputBase, parentdirs), { recursive: true });
    await fs.copyFile(path.resolve(basePath, other), path.resolve(outputBase, other));
  }
}

async function compileEmail(base, relativePath, outputBase) {
  let htmlpath = path.join(base, relativePath);
  console.log("Compiling "+htmlpath);
  let html;
  try {
    html = await fs.readFile(htmlpath);
  } catch (err) {
    if (err.code == "ENOENT") {
      console.error("File not found, skipping " + err.path);
      return;
    } else {
      throw err;
    }
  }
  let htmlstring = html.toString("utf8");
  let compiled = juice(htmlstring, {
    extraCss: await tailwindcss,
  });
  let parentdirs = path.dirname(relativePath);
  await fs.mkdir(path.resolve(outputBase, parentdirs), { recursive: true });
  await fs.writeFile(path.join(outputBase, relativePath), compiled);
}

async function listFiles(basePath) {
  let files = [];
  async function impl(basePath, relpath, result) {
    let dir = await fs.readdir(basePath, { withFileTypes: true });
    for (const file of dir) {
      if (file.isDirectory()) {
        await impl(path.join(basePath, file.name), path.join(relpath, file.name), result);
      } else if (file.isFile()) {
        result.push(path.join(relpath, file.name));
      }
    }
  }
  await impl(basePath, "", files);
  return files;
}
