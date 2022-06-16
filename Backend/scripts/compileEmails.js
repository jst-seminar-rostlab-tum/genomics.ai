const options = require("./emails.config.js");
const mjml = require("mjml");
const fs = require("fs/promises");
const path = require("path");

//const { registerComponent } = require("mjml-core");
//const MjMsoButton = require("mjml-msobutton").default;


(async function main() {
  console.log(options.input);
  const basePath = path.resolve(__dirname, options.input);
  const outputBase = path.resolve(__dirname, options.output);
  const mjmlpath = path.resolve(__dirname, options.input, "base.html");

  const basemjml = (await fs.readFile(mjmlpath)).toString("utf8");

  let dirs = await fs.readdir(basePath, { withFileTypes: true });

  console.log(dirs);
  for (const dir of dirs) {
    if(dir.isDirectory()) {
      await processFolder(basemjml, path.resolve(basePath, dir.name), path.resolve(outputBase, dir.name));
    }
  }
})();

async function processFolder(basemjml, basePath, outputBase) {
  let files = await listFiles(path.resolve(basePath));
  let others = files.filter((f) => !f.endsWith(".mjml.html"));
  console.log("Processing "+basePath);
  if(files.find((f)=>f=="content.mjml.html")) {
    await compileEmail(basemjml, basePath, outputBase);
  }
  for (const other of others) {
    console.log("Copy "+other);
    let parentdirs = path.dirname(other);
    await fs.mkdir(path.resolve(outputBase, parentdirs), { recursive: true });
    await fs.copyFile(path.resolve(basePath, other), path.resolve(outputBase, other));
  }
}

async function compileEmail(basemjml, emailpath, outputBase) {
  console.log("Compiling "+ emailpath);
  let mjmloptions = {
    filePath: emailpath,
    // mjmlConfigPath: path.join(__dirname)
  }
  let compiled = mjml(basemjml,mjmloptions)
  if(compiled.errors.length!=0) {
    console.error("Errors while processing: "+emailpath);
    compiled.errors.forEach(e=>console.error(e.formattedMessage));
  }
  await fs.mkdir(path.resolve(outputBase), { recursive: true });
  await fs.writeFile(path.join(outputBase, "index.html"), compiled.html);
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
