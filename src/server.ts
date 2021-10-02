console.log("Hello world!");

setTimeout(()=>{
    throw new Error("this is an error! did you catch me?");
}, 1500);