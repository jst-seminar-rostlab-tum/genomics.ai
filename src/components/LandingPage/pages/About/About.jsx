import React from 'react';
import styles from './about.module.css';
import { Typography } from '@mui/material';
import NavBar from '../../../NavBar/NavBar';

function About() {
  return (
    <div className={styles.headerContainer}>
      <NavBar />
      <Typography sx={{ fontWeight: '400', fontSize: '24px' }}>About us</Typography>
      <text>GeneCruncher was developed by the following team</text>
    </div>
  );
}

export default About;