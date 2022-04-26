import React, { useState } from 'react';
import styles from './tag.module.css';
import { Typography, Box } from '@mui/material';
import CancelIcon from '@mui/icons-material/Cancel';

export default function Tag(props) {
  const { text, handleClick } = props;
  const [hover, setHover] = useState(false);

  return (
    <Box
    className={styles.container}
    onMouseEnter={() => setHover(true)}
    onMouseLeave={() => setHover(false)}
    sx={{ cursor: "pointer" }}
    >
      <Typography>{text}</Typography>
      <CancelIcon
        onClick={handleClick}
        className={styles.cancelIcon}
        fontSize="smaller"
        sx={{ visibility: hover ? 'visible' : 'hidden' }}
      />
    </Box>
  );
}
