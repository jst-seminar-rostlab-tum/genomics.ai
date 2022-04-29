import * as React from 'react';
import Box from '@mui/material/Box';
import TextField from '@mui/material/TextField';
import './findMapping.css';
import SearchIcon from '@mui/icons-material/Search';

export default function FindProject() {
  return (
    // shadow
    // hight color
    // Schrift
    // color of textfield not working?
    <Box
      sx={{
        width: 350,
        justifyContent: 'center',
        marginLeft: 250,
        marginTop: 20,
        color: 'CECECE',
        borderRadius: 80,
      }}
    >
      <TextField
        className="textfield"
        id="outlined-basic"
        label={<SearchIcon />}
        variant="outlined"
        size="small"
        sx={{
          borderRadius: 100,
        }}
      />
    </Box>
  );
}
