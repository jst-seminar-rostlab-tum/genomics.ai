import React, { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import {
  Typography, createTheme, ThemeProvider, Stack, TextField,
} from '@mui/material';
import PlusIcon from 'components/GeneMapper/plusIcon';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import SearchIcon from '@mui/icons-material/Search';
import ProjectService from 'shared/services/Project.service';
import { useSubmissionProgress } from 'shared/context/submissionProgressContext';
import { MULTIPART_UPLOAD_STATUS, PROJECTS_UPDATE_INTERVAL, statusIsError } from 'shared/utils/common/constants';

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
    ProjectService.getOwnProjects().then((data) => setProjects(data));
    const timer = setInterval(() => {
      ProjectService.getOwnProjects().then((data) => setProjects(data));
      if (submissionProgress.status === MULTIPART_UPLOAD_STATUS.COMPLETE
        || submissionProgress.status === MULTIPART_UPLOAD_STATUS.CANCELING
        || statusIsError(submissionProgress.status)) {
        setSubmissionProgress({
          status: MULTIPART_UPLOAD_STATUS.IDLE,
          uploadId: '',
          chunks: 0,
          uploaded: 0,
          remaining: [],
          uploadedParts: [],
        });
      }
    }, PROJECTS_UPDATE_INTERVAL);

    return () => {
      clearInterval(timer);
    };
  }, [submissionProgress.status]);

  return (
    <div>
      <ThemeProvider theme={theme}>
        {/* {Object.entries(submissionProgress)
          .map(([key, value]) => <Typography>{`${key}: ${value}` }</Typography>)} */}
        <Box
          sx={{
            display: 'flex',
            flexDirection: 'row',
            justifyContent: 'space-between',
            pt: 2,
            pb: 3,
          }}
        >
          <Stack direction="row" className="stack" alignItems="Center">

            <Typography variant="h5" sx={{ pr: 1 }}>Your Mappings</Typography>
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
          />
        </Box>
        <div>
          {projects
            .filter((project) => (
              findString === '' || project.name.toLowerCase().includes(findString.toLowerCase())))
            .map((project) => (
              <ProjectBarCard
                key={project._id}
                projectId={project._id}
                name={project.name}
                status={project.status}
                submissionProgress={submissionProgress.uploadId === project.uploadId
                  ? submissionProgress : null}
                setSubmissionProgress={submissionProgress.uploadId === project.uploadId
                  ? setSubmissionProgress : () => {}}
              />
            ))}
        </div>

      </ThemeProvider>
    </div>
  );
}

export default GeneMapperHome;
