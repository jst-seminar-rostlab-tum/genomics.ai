import { Stack } from '@mui/material';
import React, { useCallback } from 'react';
import Header from 'components/Header';

import styles from './headerView.module.css';

function HeaderView({
  sidebarShown, title, replaceHeaderRight, children,
}) {
  const paddingL = useCallback(() => (sidebarShown ? '140px' : '350px'), [sidebarShown]);
  return (
    <Stack
      direction="column"
      sx={{
        paddingLeft: paddingL,
        paddingRight: '60px',
        width: '100%',
      }}
    >
      <Header title={title} replaceRight={replaceHeaderRight} />

      <Stack
        className="flexContainer"
        direction="row"
      >
        <div className={styles.content}>
          {children}
        </div>
      </Stack>
    </Stack>
  );
}

export default HeaderView;
