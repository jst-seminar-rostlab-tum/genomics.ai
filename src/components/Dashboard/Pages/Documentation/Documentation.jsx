import React, { useCallback } from 'react';
import { Stack, Typography } from '@mui/material';
import styles from './documentation.module.css';
import graphic1 from '../../../../assets/landing-illustrations/working.png';

function Documentation({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '130px' : '380px'), [sidebarShown]);

  return (
    <Stack
      direction="column"
      sx={{
        paddingTop: '100px',
        paddingLeft: paddingL,
        display: 'flex',
        width: '100%',
        justifyContent: 'center',
      }}
      spacing={1}
    >
      <div className={styles.title}>
        <h1>Documentation</h1>
      </div>

      <Stack direction="column" spacing={4}>
        <Typography sx={{ fontSize: '20px', width: '900px' }}>
          We are working hard on providing you with a documentation in the near future.
          If you have any problems feel free to get in touch with us.
        </Typography>
        <div>
          <img className={styles.illustration} src={graphic1} alt="working" />
        </div>
      </Stack>

    </Stack>
  );
}

export default Documentation;
