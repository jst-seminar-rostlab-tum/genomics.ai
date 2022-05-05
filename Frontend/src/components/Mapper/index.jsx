import React, { useState } from 'react';
import {
  Typography, Box, Button, IconButton, Divider, Stack, Fab,
} from '@mui/material';
import MapOutlinedIcon from '@mui/icons-material/MapOutlined';
import NoteAltIcon from '@mui/icons-material/NoteAlt';
import CancelIcon from '@mui/icons-material/Cancel';
import styles from './mapper.module.css';

function Mapper({
  initialAtlas, initialModel, editAtlasClicked, editModelClicked,
}) {
  const [atlas, setAtlas] = useState(initialAtlas);
  const [model, setModel] = useState(initialModel);
  const [open, setOpen] = useState(false);

  const handleEditAtlas = () => {
    editAtlasClicked();
    setAtlas(null);
  };

  const handleEditModel = () => {
    editModelClicked();
    setModel(null);
  };

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
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={() => setAtlas(null)}>
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
          <IconButton className={styles.iconButton} variant="contained" component="span" onClick={() => setModel(null)}>
            <CancelIcon className={styles.icon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Box className={styles.buttonBox}>
          {/* Button will be disabled if selected models and atlases are incompatible
          with each other, in this case it will be gray.
          Lets keep it enabled all the time for now. */}
          <Button className={styles.goButton} variant="contained" disableRipple>Go</Button>
        </Box>
      </Box>
      <Box className={styles.mapperBox}>
        <Fab className={styles.mapIcon} onClick={() => setOpen(!open)}>
          <MapOutlinedIcon />
        </Fab>
      </Box>
    </Box>
  );
}

export default Mapper;
