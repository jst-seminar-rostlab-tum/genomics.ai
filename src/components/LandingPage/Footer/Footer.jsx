import React from 'react';
import { Box, Typography } from '@mui/material';
import styles from './footer.module.css';
import tumLogo from '../../../assets/tum-logo.png';
import CopyrightIcon from '@mui/icons-material/Copyright';
import { fontWeight } from '@mui/system';

function Footer() {
  return (
    <Box sx={
        {padding:'25px', 
        marginTop: '250px', 
        borderTop: '1px solid black', 
        display: 'flex',
        alignItems: 'center'
        }}>
      <img className={styles.tumLogo} src={tumLogo} alt="tum logo" />
      <Typography sx={{ fontSize: '14px', padding:'10px' }}>
        <ul>
          <li><a href="about">About us</a></li>
          <li><a href="docs">Docs</a></li>
          <li><a href="contact">Contact us</a></li>
        </ul>
      </Typography>
      <CopyrightIcon sx={{paddingLeft: '20px', fontSize: '35px' }}/>
      <Typography sx={{ fontSize: '14px', padding: '10px' }}>Copyright 2021</Typography>
    </Box>
);
}

export default Footer;