import React from 'react';
import {
  Stack, Typography,
} from '@mui/material';
import styles from './docs.module.css';
import NavBar from '../../../NavBar/NavBar';
import graphic1 from '../../../../assets/landing-illustrations/working.png';
import Footer from '../../Footer/Footer';

const Docs = (props) => {
  const { setUser } = props;
  return (
    <div>
      <NavBar setUser={setUser} />

      <div className={styles.container}>
        <Stack
          spacing="30px"
        >
          <Typography sx={{ fontWeight: 'bold', fontSize: '30px' }}>Documentation</Typography>
          <div className={styles.textContainer}>
            <Typography
              sx={{ fontSize: '25px', paddingInline: '350px', maxWidth: '1700px' }}
              align="center"
            >
              We are working hard on providing you with a documentation in the near future.
              If you have any problems feel free to get in touch with us.
            </Typography>
          </div>
          <img className={styles.illustration} src={graphic1} alt="working" />
        </Stack>
      </div>

      <Footer />
    </div>

  );
};

export default Docs;
