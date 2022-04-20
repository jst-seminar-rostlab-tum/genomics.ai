import { Stack } from '@mui/material';
import React, { useCallback } from 'react';
import Header from 'components/general/Header';

import styles from './headerView.module.css';

function HeaderView({
  sidebarShown, title, rightOfTitle, replaceHeaderRight, children,
}) {
  const navbarPaddingLeft = useCallback(() => (sidebarShown ? '80px' : '190px'), [sidebarShown]);
  return (
    <Stack
      direction="column"
      sx={{
        width: '100%',
        height: '100vh',
        paddingLeft: navbarPaddingLeft,
      }}
    >
      <div style={{ paddingLeft: '60px', paddingRight: '60px' }}>
        <Header title={title} rightOfTitle={rightOfTitle} replaceRight={replaceHeaderRight} />
      </div>

      <Stack
        className="flexContainer"
        direction="row"
        sx={{ height: '100%' }}
      >
        <div className={styles.content}>
          {children}
        </div>
      </Stack>
    </Stack>
  );
}

export default HeaderView;
