import { Stack } from '@mui/material';
import React from 'react';
import Header from 'components/general/Header';

import styles from './headerView.module.css';

function HeaderView({
  title, rightOfTitle, replaceHeaderRight, children,
}) {
  return (
    <div className={styles.headerView}>
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
          sx={{ height: 'calc(100% - var(--header-height))' }}
        >
          <div className={styles.content}>
            {children}
          </div>
        </Stack>
      </Stack>
    </div>
  );
}

export default HeaderView;
