import React, { useState } from 'react';
import Box from '@mui/material/Box';
import FindMapping from 'components/GeneMapper/findMapping';
import {
  Typography, createTheme, ThemeProvider, Stack,
} from '@mui/material';
import PlusIcon from 'components/GeneMapper/plusIcon';
import styles from './geneMapper.css';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import FindMappingRequest from 'components/GeneMapper/findMappingRequest';

const theme = createTheme({
  palette: {
    primary: {
      main: '#000000',
    },
    secondary: {
      main: '#5676E4',
    },
  },
  typography: {
    fontFamily: 'Lato',
    fontSize: 16,
    // the size of save
  },
});
theme.typography.h3 = {
  fontSize: '1.2rem',
  '@media (min-width:600px)': {
    fontSize: '1.5rem',
  },
  [theme.breakpoints.up('md')]: {
    fontSize: '2rem',
  },
};
const themeIcon = createTheme({
  palette: {
    primary: {
      main: '#5676E4',
    },
  },
});

function GeneMapperHome() {
  return (
    <div>
      <ThemeProvider theme={theme}>
        <Box className="box">
          <Stack direction="row" spacing={160} sx={{ marginTop: 15, marginLeft: 30 }} className="stack">
            <Stack direction="row" spacing={4} className="stack">
              <Typography variant="h3" color="primary" className={styles.title}>Your Mappings </Typography>
              <ThemeProvider theme={themeIcon}>
                <PlusIcon />
              </ThemeProvider>
            </Stack>
            <FindMappingRequest className="findMappingRequest" />
          </Stack>
          <ProjectBarCard />
          <ProjectBarCard />
          <ProjectBarCard />
        </Box>
      </ThemeProvider>

    </div>
  );
}

export default GeneMapperHome;
