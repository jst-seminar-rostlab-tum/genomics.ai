import { buttonWidth, cornerRadius, fill } from "./constants";

const addSmallButton = (container) => {
  container
    .append("rect")
    .attr("fill", "lightgrey")
    .attr("width", 15)
    .attr("height", 15)
    .attr("stroke", "black")
    .attr("stroke-width", 2)
    .attr("rx", cornerRadius);

  container.on("mouseover", () => {
    container.select("rect").transition().attr("fill", "whitesmoke");
  });

  container.on("mouseleave", () => {
    container.select("rect").transition().attr("fill", "lightgrey");
  });
};

export const addSmallCloseButton = (container) => {
  addSmallButton(container);

  container
    .append("line")
    .attr("x1", 3)
    .attr("x2", 12)
    .attr("y1", 3)
    .attr("y2", 12)
    .attr("stroke", "black")
    .attr("stroke-width", 1);

  container
    .append("line")
    .attr("x1", 12)
    .attr("x2", 3)
    .attr("y1", 3)
    .attr("y2", 12)
    .attr("stroke", "black")
    .attr("stroke-width", 1);
};

const addButton = (container, height = buttonWidth) => {
  const rect = container
    .append("rect")
    .attr("width", buttonWidth)
    .attr("height", height)
    .attr("fill", fill)
    .attr("rx", cornerRadius)
    .attr("stroke", fill)
    .attr("stroke-width", 1);

  container.on("mouseover", () => {
    rect.transition().attr("fill", "dimgray");
  });

  container.on("mouseleave", () => {
    rect.transition().attr("fill", fill);
  });

  container.on("click.format", () => {
    rect.transition().attr("fill", "grey").transition().attr("fill", "dimgray");
  });
};

const addZoomInButton = (container) => {
  addButton(container);

  container
    .append("line")
    .attr("x1", buttonWidth / 2)
    .attr("y1", 8)
    .attr("x2", buttonWidth / 2)
    .attr("y2", buttonWidth - 8)
    .attr("stroke", "lightgrey")
    .attr("stroke-width", 4);

  container
    .append("line")
    .attr("x1", 8)
    .attr("y1", buttonWidth / 2)
    .attr("x2", buttonWidth - 8)
    .attr("y2", buttonWidth / 2)
    .attr("stroke", "lightgrey")
    .attr("stroke-width", 4);
};

const addZoomOutButton = (container) => {
  addButton(container);

  container
    .append("line")
    .attr("x1", 8)
    .attr("y1", buttonWidth / 2)
    .attr("x2", buttonWidth - 8)
    .attr("y2", buttonWidth / 2)
    .attr("stroke", "lightgrey")
    .attr("stroke-width", 4);
};

const addToggleLabelsButton = (container) => {
  addButton(container, buttonWidth / 2);

  container
    .append("text")
    .attr("fill", "lightgrey")
    .attr("font-size", 11)
    .attr("x", 4)
    .attr("y", 14)
    .text("Labels");
};

const addResetButton = (container) => {
  addButton(container, buttonWidth / 2);

  container
    .append("text")
    .attr("fill", "lightgrey")
    .attr("font-size", 11)
    .attr("x", 6)
    .attr("y", 14)
    .text("Reset");
};

export const addControlButtons = (container) => {
  const zoomIn = container.append("g");
  const zoomOut = container
    .append("g")
    .attr("transform", `translate(0,${buttonWidth + 5})`);
  const toggleLabels = container
    .append("g")
    .attr("transform", `translate(0,${2 * buttonWidth + 2 * 5})`);
  const reset = container
    .append("g")
    .attr("transform", `translate(0,${2.5 * buttonWidth + 3 * 5})`);

  addZoomInButton(zoomIn);
  addZoomOutButton(zoomOut);
  addToggleLabelsButton(toggleLabels);
  addResetButton(reset);

  return { zoomIn, zoomOut, toggleLabels, reset };
};
