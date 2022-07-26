import React, { useState, useEffect, useContext } from 'react';
import {
  Typography, Box, IconButton, Divider, Stack, Fab,
} from '@mui/material';
import MapOutlinedIcon from '@mui/icons-material/MapOutlined';
import CancelIcon from '@mui/icons-material/Cancel';
import styles from './mapper.module.css';
import CustomButton from 'components/CustomButton';
import { LoginContext } from 'shared/context/loginContext';
import { FaExclamationTriangle } from "react-icons/fa";

function Mapper({
  mapperAtlas, mapperModel, handleAtlasSelection,
  handleModelSelection, open, fabOnClick, handleMap,
}) {
  const [atlas, setAtlas] = useState(mapperAtlas);
  const [model, setModel] = useState(mapperModel);

  const deleteAtlas = () => {
    handleAtlasSelection(null);
    setAtlas(null);
  };

  const deleteModel = () => {
    handleModelSelection(null);
    setModel(null);
  };

  useEffect(() => {
    setAtlas(mapperAtlas);
    setModel(mapperModel);
  }, [mapperAtlas, mapperModel]);

  return (
    <Box className={styles.container}>
      <Box className={styles.borderContainer} sx={{ display: open ? 'grid' : 'none', boxShadow: '0px 4px 6px rgba(0, 0, 0, 0.25)' }}>
        <Typography className={styles.title}>Mapper</Typography>
        <Divider className={styles.divider} />
        <Typography>Selected Atlas</Typography>
        <Stack
          direction="row"
        >
          <Typography className={styles.filename}>{atlas || 'None'}</Typography>
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={deleteAtlas} disabled={!atlas}>
            <CancelIcon className={atlas ? styles.icon : styles.disabledIcon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Typography>Selected Model</Typography>
        <Stack
          direction="row"
        >
          <Typography className={styles.filename}>{model || 'None'}</Typography>
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={deleteModel} disabled={!model}>
            <CancelIcon className={model ? styles.icon : styles.disabledIcon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Box className={styles.buttonBox}>
          <CustomButton type="primary" disabled={!model || !atlas} onClick={handleMap}>Create Project</CustomButton>
        </Box>
      </Box>
      <Box className={styles.mapperBox}>
        <Fab className={styles.mapIcon} onClick={fabOnClick}>
          <MapOutlinedIcon />
        </Fab>
      </Box>
    </Box>
  );
}

export default Mapper;
