import React, { useState, useEffect } from 'react';
import {
  Typography, Box, IconButton, Divider, Stack, Fab,
} from '@mui/material';
import MapOutlinedIcon from '@mui/icons-material/MapOutlined';
import CancelIcon from '@mui/icons-material/Cancel';
import styles from './mapper.module.css';
import CustomButton from 'components/CustomButton';
import {
  useHistory, useRouteMatch,
} from 'react-router-dom';

function Mapper({
  mapperAtlas, mapperModel, handleAtlasSelection, handleModelSelection, open, fabOnClick,
  fabVisible,
}) {
  const [atlas, setAtlas] = useState(mapperAtlas);
  const [model, setModel] = useState(mapperModel);
  const { path } = useRouteMatch();
  const history = useHistory();

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
    <Box className={styles.container} sx={{ display: fabVisible ? 'block' : 'none' }}>
      <Box className={styles.borderContainer} sx={{ display: open ? 'block' : 'none', boxShadow: '0px 4px 6px rgba(0, 0, 0, 0.25)' }}>
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
          {/* Button will be disabled if selected models and atlases are incompatible with eachother, in this case it will be gray. Lets keep it enabled all the time for now. */}
          <CustomButton disabled={!atlas || !model} type="primary" onClick={() => { history.push(`${path}/genemapper/create`); }}>Go</CustomButton>
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
