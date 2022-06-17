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
        minHeight: "calc(100vh - 40px)",
        padding: "12px 60px 0 60px", 
        boxSizing: "border-box",
        display: "flex",
        }}>
        <Stack
          direction="column"
          sx={{
            // height: "98vh",
            width: '100%',
            paddingLeft: '80px',
            flexGrow: 1
          }}
        >
          <div style={{ paddingLeft: '60px', paddingRight: '60px' }}>
            <Header title={title} rightOfTitle={rightOfTitle} replaceRight={replaceHeaderRight} />
          </div>

          <Stack ref={ref}
            className="flexContainer"
            direction="row"
            sx={{ flexGrow: 1 }}
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
