import { Box } from "@mui/system";
import csv from "assets/test_file2.csv";
import ResultVisualization from "components/GeneMapper/ResultVisualization";
import { useParams } from "react-router-dom";

export default function() {

  const { id } = useParams()

  return (
    <Box>
      <ResultVisualization dataUrl={csv} onlyUmap/>
    </Box>
  )
}