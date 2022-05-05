import React, { useEffect, useState } from 'react';
import { Box, Typography } from '@mui/material';
import CustomButton from 'components/CustomButton';
import AtlasService from 'shared/services/Atlas.service';
import { useParams, useHistory, useLocation } from 'react-router-dom';

export const LearnMoreAtlasComponent = ({ onClick, id }) => {
  const [atlas, setAtlas] = useState(null);

  useEffect(() => {
    if (id) {
      AtlasService.getAtlasById(id)
        .then((data) => setAtlas(data))
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
        <Typography sx={{ fontSize: '36px', fontWeigth: 700 }}>{atlas?.name}</Typography>
      </Box>
      <Box>
        <Typography sx={{ fontSize: '20px', fontWeight: 600, borderBottom: '1px solid black' }}>Overview</Typography>
      </Box>
      <Box sx={{ display: 'flex', flexDirection: 'row', paddingTop: '16px' }}>
        <Typography sx={{ fontSize: '16px', fontWeight: 500 }}>
          Modalities:
          &nbsp;
        </Typography>
        <Typography sx={{ fontSize: '16px', fontWeight: 300 }}>
          {atlas?.modalities}
        </Typography>
      </Box>
      <Box sx={{ display: 'flex', flexDirection: 'row' }}>
        <Typography sx={{ fontSize: '16px', fontWeight: 500 }}>
          Cells in Reference:
          &nbsp;
        </Typography>
        <Typography sx={{ fontSize: '16px', fontWeight: 300 }}>
          {atlas?.numberOfCells}
        </Typography>
      </Box>
      <Box sx={{ display: 'flex', flexDirection: 'row' }}>
        <Typography sx={{ fontSize: '16px', fontWeight: 500 }}>
          Species:
          &nbsp;
        </Typography>
        <Typography sx={{ fontSize: '16px', fontWeight: 300 }}>
          {atlas?.species}
        </Typography>
      </Box>

      <CustomButton sx={{ marginTop: "1em", padding: "0.5em 2em 0.5em 2em" }} type="primary" onClick={onClick}>Map</CustomButton>
    </Box>
  );
};

export default function LearnMore() {
  const history = useHistory();
  const path = useLocation();
  const { id } = useParams();
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
      <Box sx={{
        marginTop: '12px',
        width: '80%',
        minWidth: '1200px',
      }}
      >
        <LearnMoreAtlasComponent id={id} onClick={() => history.push(`${path.pathname}/visualization`)} />
      </Box>
    </Box>
  );
}
