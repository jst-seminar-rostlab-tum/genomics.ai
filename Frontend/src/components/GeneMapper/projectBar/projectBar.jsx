import * as React from 'react';
import Box from '@mui/material/Box';
import Button from '@mui/material/Button';

export default function ProjectBar() {
  return (
    <div>

      <Box
        component="span"
        sx={{
          p: 5, border: '1px dashed grey', width: 350, length: 1000,
        }}
      >
        <Button>Save</Button>
      </Box>
    </div>

  );
}
/*
sx={{
    '& > :not(style)': { m: 1, width: '25ch' },
    width: 350,
    justifyContent: 'center',
    position: 'absolute',
    top: '20%',
    right: '5%',
    color: 'CECECE',
  }}
  */
