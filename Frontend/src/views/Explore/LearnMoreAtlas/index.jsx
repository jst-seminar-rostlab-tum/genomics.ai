import React, { useEffect, useState } from 'react';
import { Box, Typography } from '@mui/material';
import CustomButton from 'components/CustomButton';
import AtlasService from 'shared/services/Atlas.service';
import { useParams, useHistory, useLocation } from 'react-router-dom';
import { colors } from 'shared/theme/colors';

export const LearnMoreAtlasComponent = ({ onClick, id, isMap = false, isSearchPage = false }) => {
  const [atlas, setAtlas] = useState(null);
  const history = useHistory();
  const path = useLocation().pathname;

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
      {isSearchPage && <Typography onClick={history.goBack} sx={{ cursor: "pointer", fontSize: "18px", fontWeight: 500,  color: colors.neutral[800], ":hover": { color: colors.primary[500] }}}>Go Back</Typography>}
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
      {
        isMap
        &&
        !isSearchPage 
        &&
        <>
          <CustomButton sx={{ marginTop: '1em', padding: '0.5em 2em 0.5em 2em' }} type="secondary" onClick={() => history.push(`${path}/visualization`)}>Visualize</CustomButton>
          <CustomButton sx={{ marginTop: '1em', padding: "0.5em 2em 0.5em 2em" }} type="primary" onClick={() => onClick(atlas)}>Select</CustomButton>
        </>
      }
    </Box>
  );
};

export default function LearnMore({ handleSelect }) {
  const path = useLocation();
  const { id } = useParams();
  
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
      <Box sx={{
        marginTop: '12px',
        width: '80%',
      }}
      >
        <LearnMoreAtlasComponent id={id} isMap={true} onClick={handleSelect} isSearchPage={path.pathname.includes("search")}/>
      </Box>
    </Box>
  );
}
