import React from 'react';
import { Box, Typography } from '@mui/material';
import CopyrightIcon from '@mui/icons-material/Copyright';
import { Link } from 'react-router-dom';
import styles from './footer.module.css';
import tumLogo from '../../../assets/tum-logo.png';

function Footer() {
  return (
    <Box sx={
        {
          padding: '25px',
          marginTop: '250px',
          borderTop: '1px solid black',
          display: 'flex',
          alignItems: 'center',
        }
}
    >
      <img className={styles.tumLogo} src={tumLogo} alt="tum logo" />
      <Typography sx={{ fontSize: '14px', padding: '10px' }}>
        <ul>
          <li><Link to="about">About us</Link></li>
          <li><Link to="docs">Docs</Link></li>
          <li><Link to="contact">Contact us</Link></li>
        </ul>
      </Typography>
      <CopyrightIcon sx={{ paddingLeft: '20px', fontSize: '35px' }} />
      <Typography sx={{ fontSize: '14px', padding: '10px' }}>Copyright 2021</Typography>
    </Box>
  );
}

export default Footer;
