import React, { useState } from 'react';
import styles from './tag.module.css';
import { Typography, Box } from '@mui/material';
import CancelIcon from '@mui/icons-material/Cancel';
import { borderRadius } from '@mui/system';

export default function Tag(props) {
  const { text, handleClick } = props;
  const [hover, setHover] = useState(false);

  return (
    <Box
      onMouseEnter={() => setHover(true)}
      onMouseLeave={() => setHover(false)}
      onClick={handleClick}
      sx={{
        backgroundColor: '#5576E4', borderRadius: 4, paddingX: 2, paddingY: 1, color: 'white', marginX: 0.5, position: 'relative',
      }}
    >
      <Typography>{text}</Typography>
      <CancelIcon sx={{ visibility: hover ? 'visible' : 'hidden', color: 'red', fontSize: 'small', position: 'absolute', zIndex: 2, top: 0, right: 0 }} />
    </Box>
  );
}
