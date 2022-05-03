import { Box } from "@mui/system";
import ResultVisualization from "components/GeneMapper/ResultVisualization";
import { useParams } from "react-router-dom";

export default function() {

  const { id } = useParams

  return (
    <Box>
      <ResultVisualization dataUrl={path} onlyUmap/>
    </Box>
  )
}