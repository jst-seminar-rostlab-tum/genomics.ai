import * as React from 'react';
import Box from '@mui/material/Box';
import TextField from '@mui/material/TextField';
import './findMapping.css';
import SearchIcon from '@mui/icons-material/Search';
import { Stack, Typography } from '@mui/material';

export default function FindProject() {
  return (
    // shadow
    // hight color
    // Schrift
    // color of textfield not working?
    <Box>
      <TextField
        className="textfield"
        id="outlined-basic"
        label={(
          <Stack>
            <SearchIcon />
            Find a Mapping
          </Stack>
        )}
        variant="outlined"
        size="small"
        value={findString}
      />
    </Box>
  );
}
