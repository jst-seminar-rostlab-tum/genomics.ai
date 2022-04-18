import React, { useState } from 'react';
import {
  Typography, Box, Button, Input, IconButton, Divider, Stack,
} from '@mui/material';
import MapIcon from '@mui/icons-material/Map';
import NoteAltIcon from '@mui/icons-material/NoteAlt';
import CancelIcon from '@mui/icons-material/Cancel';

function Mapper() {
  const [files, setFiles] = useState({
    atlas: null,
    model: null,
  });

  return (
    <Box>
      <Typography>Mapper</Typography>
      <Divider />
      <Typography>Selected Atlas</Typography>
      <Stack
        direction="row"
      >
        <Typography>{files.atlas ? files.atlas.name : 'None'}</Typography>
        <label htmlFor="select-atlas">
          <Input id="select-atlas" type="file" sx={{ display: 'none' }} onChange={(e) => setFiles({ ...files, atlas: e.target.files[0] })} />
          <IconButton variant="contained" component="span">
            <NoteAltIcon />
          </IconButton>
        </label>
        <IconButton variant="contained" component="span" onClick={() => setFiles({ ...files, atlas: null })}>
          <CancelIcon />
        </IconButton>
      </Stack>
      <Divider />
      <Typography>Selected Model</Typography>
      <Stack
        direction="row"
      >
        <Typography>{files.atlas ? files.model.name : 'None'}</Typography>
        <label htmlFor="select-model">
          <Input id="select-model" type="file" sx={{ display: 'none' }} onChange={(e) => setFiles({ ...files, model: e.target.files[0] })} />
          <IconButton variant="contained" component="span">
            <NoteAltIcon />
          </IconButton>
        </label>
        <IconButton variant="contained" component="span" onClick={() => setFiles({ ...files, model: null })}>
          <CancelIcon />
        </IconButton>
      </Stack>
      <Divider />
      <Button variant="contained">Go</Button>
      <MapIcon />
    </Box>
  );
}

export default Mapper;
