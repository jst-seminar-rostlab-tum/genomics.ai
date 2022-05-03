import React, { useEffect, useState } from 'react';
import { Box, Typography } from '@mui/material';
import ModelsService from 'shared/services/Models.service';
import { useParams } from 'react-router-dom';

export const LearnMoreModelComponent = () => {
  const { id } = useParams()
  const [model, setModel] = useState(null);
  
  useEffect(() => {
    if(id) {
      ModelsService.getModelById(id)
        .then((data) => setModel(data))
        .catch((err) => console.log(err));
    }
  }, [id]);
  
  return (
    <Box sx={{
      display: 'flex',
      flexDirection: 'column',
      alignItems: 'flex-start',
      justifyContent: 'space-between',
    }}
    >
      <Box sx={{
        display: 'flex', flexDirection: 'row', width: '100%', justifyContent: 'space-between',
      }}
      >
        <Typography sx={{ fontSize: '36px', fontWeigth: 700 }}>{model?.name}</Typography>
      </Box>
      <Box>
        <Typography sx={{ fontSize: '20px', fontWeight: 600, borderBottom: '1px solid black' }}>Overview</Typography>
      </Box>
      
      <Typography>Description: {model?.description}</Typography> 
    </Box>
  );
}

export default function LearnMore() {

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
      <Box sx={{
        marginTop: '12px',
        width: '80%',
        minWidth: '1200px',
      }}
      >
        <LearnMoreModelComponent />
      </Box>
    </Box>
  )
}
