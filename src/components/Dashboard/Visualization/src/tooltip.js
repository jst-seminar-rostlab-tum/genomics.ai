import * as d3 from "d3";

const opacity = 0.8;

export const addTooltip = (tooltipDiv, dots, clusterMapping) => {
  const showTooltip = function (d) {
    tooltipDiv.style("display", "block")
    d3.select(this).style("stroke", "black")
  };

  const setToolTipText = (e, d) => {
    const info = `Id: ${d.id}\nCell-Type: ${clusterMapping.get(d.cluster)}`;
    tooltipDiv.style("width", `${28 + d.id.length * 4}pt`)
    tooltipDiv.html(info)
    .style("position", "absolute")
    .style("left", `${e.clientX}px`)
    .style("top", `${e.clientY}px`);
  };

  const hideToolTip = function (d) {
    tooltipDiv.style("display", "none");
    d3.select(this)
        .style("stroke", "none")
        .style("opacity", opacity);
  };

  dots
    .on("mouseover", showTooltip)
    .on("mousemove", setToolTipText)
    .on("mouseleave", hideToolTip);

  return tooltipDiv;
}