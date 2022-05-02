import React, { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import {
  Typography, createTheme, ThemeProvider, Stack, TextField,
} from '@mui/material';
import PlusIcon from 'components/GeneMapper/plusIcon';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import SearchIcon from '@mui/icons-material/Search';
import ProjectMock from 'shared/services/mock/projects';
import { useSubmissionProgress } from 'shared/context/submissionProgressContext';
import { getSubmissionProgressPercentage } from 'shared/services/UploadLogic';

const theme = createTheme({
  palette: {
    primary: {
      main: '#000000',
    },
    secondary: {
      main: '#5676E4',
    },
  },
});

const themeIcon = createTheme({
  palette: {
    primary: {
      main: '#5676E4',
    },
  },
});

function GeneMapperHome() {
  const [projects, setProjects] = useState([]);
  const [findString, setFindString] = useState('');
  const [submissionProgress, setSubmissionProgress] = useSubmissionProgress();

  useEffect(() => {
    console.log(getSubmissionProgressPercentage(submissionProgress));
  }, [submissionProgress]);

  useEffect(() => {
    ProjectMock.getProjects().then((data) => setProjects(data));
  }, []);

  return (
    <div>
      <ThemeProvider theme={theme}>
        <Box
          sx={{
            display: 'flex',
            flexDirection: 'row',
            justifyContent: 'space-between',
            paddingBottom: '1em',
            marginTop: '1px',
          }}
        >
          <Stack direction="row" className="stack">

            <Box sx={{ fontSize: 'h6.fontSize', padding: 3 }}>Your Mappings</Box>
            <ThemeProvider theme={themeIcon}>
              <PlusIcon />
            </ThemeProvider>
          </Stack>
          <TextField
            id="outlined-basic"
            label={(
              <Stack direction="row">
                <SearchIcon />
                Find a Mapping
              </Stack>
        )}
            variant="outlined"
            size="small"
            value={findString}
            onChange={(e) => setFindString(e.target.value)}
            sx={{ m: 2.5 }}
          />
        </Box>
        <div>
          {projects
            .filter((project) => (
              findString === '' || project.name.toLowerCase().includes(findString.toLowerCase())))
            .map((project) => (
              <ProjectBarCard projectId={project._id} name={project.name} status={project.status} />
            ))}
        </div>

      </ThemeProvider>
    </div>
  );
}

export default GeneMapperHome;
