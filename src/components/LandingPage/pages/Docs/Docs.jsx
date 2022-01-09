import React from 'react';
import {
  Stack, Typography,
} from '@mui/material';
import styles from './docs.module.css';
import NavBar from '../../../NavBar/NavBar';
import graphic1 from '../../../../assets/landing-illustrations/working.png';

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
          <Typography sx={{ fontSize: '25px', paddingInline: '350px' }}>
            We are working hard on providing you with a documentation in the near future.
            If you have any problems feel free to get in touch with us.
          </Typography>
          <img className={styles.illustration} src={graphic1} alt="working" />
        </Stack>
      </div>
    </div>

  );
};

export default Docs;
