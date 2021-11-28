import React from 'react';
import NavBar from '../../../NavBar/NavBar';
import styles from './contact.module.css';
import { TextField} from '@mui/material';

function Contact() {
  return (
      <div className={styles.headerContainer}>
          <NavBar />
          <h1>Contact us</h1>
          
      </div>
  )
}

export default Contact;