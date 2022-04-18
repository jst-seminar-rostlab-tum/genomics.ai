import React, { useState } from 'react';
import {
  Typography, Box, Button, Input, IconButton, Divider, Stack,
} from '@mui/material';
import MapOutlinedIcon from '@mui/icons-material/MapOutlined';
import NoteAltIcon from '@mui/icons-material/NoteAlt';
import CancelIcon from '@mui/icons-material/Cancel';
import styles from './mapper.module.css';

function Mapper() {
  const [files, setFiles] = useState({
    atlas: null,
    model: null,
  });

  return (
    <Box className={styles.container}>
      <Box className={styles.borderContainer} sx={{ boxShadow: 5 }}>
        <Typography className={styles.title}>Mapper</Typography>
        <Divider className={styles.divider} />
        <Typography>Selected Atlas</Typography>
        <Stack
          direction="row"
        >
          <Typography className={styles.filename}>{files.atlas ? files.atlas.name : 'None'}</Typography>
          <label htmlFor="select-atlas">
            <Input id="select-atlas" type="file" className={styles.input} onChange={(e) => setFiles({ ...files, atlas: e.target.files[0] })} />
            <IconButton variant="contained" component="span">
              <NoteAltIcon className={styles.icon} />
            </IconButton>
          </label>
          <IconButton variant="contained" component="span" onClick={() => setFiles({ ...files, atlas: null })}>
            <CancelIcon className={styles.icon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Typography>Selected Model</Typography>
        <Stack
          direction="row"
        >
          <Typography className={styles.filename}>{files.model ? files.model.name : 'None'}</Typography>
          <label htmlFor="select-model">
            <Input id="select-model" type="file" className={styles.input} onChange={(e) => setFiles({ ...files, model: e.target.files[0] })} />
            <IconButton variant="contained" component="span">
              <NoteAltIcon className={styles.icon} />
            </IconButton>
          </label>
          <IconButton variant="contained" component="span" onClick={() => setFiles({ ...files, model: null })}>
            <CancelIcon className={styles.icon} />
          </IconButton>
        </Stack>
        <Divider className={styles.divider} />
        <Box className={styles.buttonBox}>
          <Button className={styles.goButton} variant="contained">Go</Button>
        </Box>
      </Box>
      <Box className={styles.mapperBox}>
        <IconButton className={styles.mapIcon}>
          <MapOutlinedIcon />
        </IconButton>
      </Box>
    </Box>
  );
}

export default Mapper;
