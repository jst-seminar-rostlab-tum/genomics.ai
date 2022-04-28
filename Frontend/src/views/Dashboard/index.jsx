import { Container, Stack } from '@mui/material';
import React, { useCallback } from 'react';
import Uploader from 'components/Uploader';
import StatusQueue from 'components/StatusQueue';
import styles from './dashboard.module.css';

import { TabGroup } from 'components/Tab'
import {ModelCard} from 'components/Cards/ModelCard'
import { useState } from 'react'
import {Box} from '@mui/material'

const x=(<Box sx={{position: "relative", padding: "10px", display: "flex", flexDirection: "column", justifyContent: "left", gap: "5px"}}>
<ModelCard title="model 1" description="asdfasdfasdfasdfasdf"/>
<ModelCard title="model 1" description="asdfasdfasdfasdfasdf"/>
<ModelCard title="model 1" description="asdfasdfasdfasdfasdf"/>
<ModelCard title="model 1" description="asdfasdfasdfasdfasdf"/></Box>
)

function Dashboard({ sidebarShown }) {
  const [value, setValue] = useState(0)

  const paddingL = useCallback(() => (sidebarShown ? '100px' : '350px'), [sidebarShown]);
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
    >
      <div className={styles.title}>
        <h1>Dashboard</h1>
      </div>

      <Stack
        className="flexContainer"
        direction="row"
      >
        <Container className={styles.fileUpload}>
          <Uploader />
        </Container>
        <Container className={styles.fileQueue}>
          <StatusQueue />
        </Container>
      </Stack>
      
      <TabGroup value={value} setValue={setValue} width="400px" height="400px" tabsInfo={[{label: "tab1", additionalContent: x},{label: "tab2"},{label: "tab3"}]} />
    </Stack>
  );
}

export default Dashboard;
