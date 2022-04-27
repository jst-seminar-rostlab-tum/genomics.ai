import React, { useState, useEffect } from 'react';
import {
  Typography, Box, Button, IconButton, Divider, Stack, Fab,
} from '@mui/material';
import MapOutlinedIcon from '@mui/icons-material/MapOutlined';
import NoteAltIcon from '@mui/icons-material/NoteAlt';
import CancelIcon from '@mui/icons-material/Cancel';
import styles from './mapper.module.css';
import { useHistory } from 'react-router-dom';

function Mapper({
  mapperAtlas, mapperModel, setSelectedAtlas, setSelectedModel, open, fabOnClick,
}) {
  const [atlas, setAtlas] = useState(mapperAtlas);
  const [model, setModel] = useState(mapperModel);
  const history = useHistory();

  const handleEditAtlas = () => {
    history.push('/explore/atlases');
    setSelectedAtlas(-1);
    setAtlas(null);
  };

  const deleteAtlas = () => {
    setSelectedAtlas(-1);
    setAtlas(null);
  };

  const handleEditModel = () => {
    history.push('/explore/models');
    setSelectedModel(-1);
    setModel(null);
  };

  const deleteModel = () => {
    setSelectedModel(-1);
    setModel(null);
  };

  useEffect(() => {
    setAtlas(mapperAtlas);
    setModel(mapperModel);
  }, [mapperAtlas, mapperModel]);

  return (
    <Box className={styles.container}>
      <Box className={styles.borderContainer} sx={{ visibility: open ? 'visible' : 'hidden', boxShadow: '0px 4px 6px rgba(0, 0, 0, 0.25)' }}>
        <Typography className={styles.title}>Mapper</Typography>
        <Divider className={styles.divider} />
        <Typography>Selected Atlas</Typography>
        <Stack
          direction="row"
        >
          <Typography className={styles.filename}>{atlas || 'None'}</Typography>
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={handleEditAtlas}>
            <NoteAltIcon className={styles.icon} />
          </IconButton>
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={deleteAtlas}>
            <CancelIcon className={styles.icon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Typography>Selected Model</Typography>
        <Stack
          direction="row"
        >
          <Typography className={styles.filename}>{model || 'None'}</Typography>
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={handleEditModel}>
            <NoteAltIcon className={styles.icon} />
          </IconButton>
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={deleteModel}>
            <CancelIcon className={styles.icon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Box className={styles.buttonBox}>
          {/* Button will be disabled if selected models and atlases are incompatible with eachother, in this case it will be gray. Lets keep it enabled all the time for now. */}
          <Button className={styles.goButton} variant="contained" disableRipple>Go</Button>
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
