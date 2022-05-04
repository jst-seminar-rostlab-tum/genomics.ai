import { CircularProgress } from '@mui/material';
import { Box } from '@mui/system';
import ResultVisualization from 'components/GeneMapper/ResultVisualization';
import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import axiosInstance from 'shared/services/axiosInstance';

export default function () {
  const { id } = useParams();
  const [data, setData] = useState(null);

  useEffect(() => {
    axiosInstance
      .get(`/atlas/${id}/visualization`)
      .then((data) => setData(data));
  }, []);
  console.log(data);
  return (
    <Box sx={{ height: 500 }}>
      {!data && <CircularProgress />}
      {data && <ResultVisualization dataUrl={data.data} onlyUmap />}
    </Box>
  );
}
