import { Box, Stack } from '@mui/material';
import React, { useEffect, useRef } from 'react';
import Header from 'components/general/Header';
import Footer from 'components/Footer';

import styles from './headerView.module.css';

function HeaderView({
  title, rightOfTitle, replaceHeaderRight, children,
}) {

  const ref=useRef();

  useEffect(()=>{
    console.log(ref)
  }, [])

  return (
    <Box ref={ref} sx={{
      width: "100%", 
      height: "100vh",
      overflowY: "auto"}}>
    {/* <div className={styles.headerView}> */}
      <Box sx={{
        minHeight: "92vh",
        padding: "12px 60px 0 60px", 
        boxSizing: "border-box"}}>
        <Stack
          direction="column"
          sx={{
            // height: "98vh",
            width: '100%',
            paddingLeft: '80px',
          }}
        >
          <div style={{ paddingLeft: '60px', paddingRight: '60px' }}>
            <Header title={title} rightOfTitle={rightOfTitle} replaceRight={replaceHeaderRight} />
          </div>

          <Stack ref={ref}
            className="flexContainer"
            direction="row"
            sx={{ height: 'calc(92vh - var(--header-height))' }}
          >
            <div className={styles.content}>
              {children}
            </div>
          </Stack>
        </Stack>
      </Box>
      <Footer />
    </Box>
    // </div>
  );
}

export default HeaderView;
