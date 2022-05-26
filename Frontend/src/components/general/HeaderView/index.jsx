import { Box, Stack } from '@mui/material';
import React from 'react';
import Header from 'components/general/Header';

import styles from './headerView.module.css';

function HeaderView({
  title, rightOfTitle, replaceHeaderRight, children,
}) {
  return (
    // <div className={styles.headerView}>
    <Box sx={{
      width: "100%", 
      height: "92vh", 
      maxHeight: "92vh", 
      overflowY: "auto", 
      padding: "12px 60px 0 60px", 
      boxSizing: "border-box"}}>
      <Stack
        direction="column"
        sx={{
          width: '100%',
          height: '100vh',
          paddingLeft: '80px',
        }}
      >
        <div style={{ paddingLeft: '60px', paddingRight: '60px' }}>
          <Header title={title} rightOfTitle={rightOfTitle} replaceRight={replaceHeaderRight} />
        </div>

        <Stack
          className="flexContainer"
          direction="row"
          // sx={{ height: 'calc(92vh - var(--header-height))' }}
        >
          <div className={styles.content}>
            {children}
          </div>
        </Stack>
      </Stack>
    </Box>
    // </div>
  );
}

export default HeaderView;
