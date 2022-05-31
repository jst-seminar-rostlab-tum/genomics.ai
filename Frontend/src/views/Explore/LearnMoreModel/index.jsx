import React, { useEffect, useState } from 'react';
import { Box, Link, Typography } from '@mui/material';
import ModelsService from 'shared/services/Models.service';
import { useParams } from 'react-router-dom';
import CustomButton from 'components/CustomButton';
import { useHistory, useLocation } from 'react-router-dom/cjs/react-router-dom.min';
import { colors } from "shared/theme/colors";

export const LearnMoreModelComponent = ({ onClick, id, isMap = false, isSearchPage = false }) => {

  const [model, setModel] = useState(null);
  const history = useHistory();

  useEffect(() => {
    if (id) {
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
      {isSearchPage && <Typography onClick={history.goBack} sx={{ cursor: "pointer", fontSize: "18px", fontWeight: 500,  color: colors.neutral[800], ":hover": { color: colors.primary[500] }}}>Go Back</Typography>}
      <Box sx={{
        display: 'flex', flexDirection: 'row', width: '100%', justifyContent: 'space-between',
      }}
      >
        <Typography sx={{ fontSize: '36px', fontWeigth: 700 }}>{model?.name}</Typography>
      </Box>
      <Box>
        <Typography sx={{ fontSize: '20px', fontWeight: 600, borderBottom: '1px solid black' }}>Overview</Typography>
      </Box>

      <Typography sx={{ width: '100%', maxWidth: '800px' }}>Description: {model?.description}</Typography>
      {
        // isSelect
        isMap && !isSearchPage &&
        <CustomButton sx={{ marginTop: "1em", padding: "0.5em 2em 0.5em 2em" }} type="primary" onClick={() => onClick(model)}>Select</CustomButton>
      }
    </Box>
  );
}

export default function LearnMore({handleSelect}) {
  const { id } = useParams();
  const path = useLocation();

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
      <Box sx={{
        marginTop: '12px',
        width: '80%',
      }}
      >
        <LearnMoreModelComponent id={id} isMap={true} onClick={handleSelect} isSearchPage={path.pathname.includes("search")}/>
      </Box>
    </Box>
  )
}
