import React, { useState } from 'react';
import Box from '@mui/material/Box';
import {
  Typography, createTheme, ThemeProvider, Stack, TextField,
} from '@mui/material';
import PlusIcon from 'components/GeneMapper/plusIcon';
import styles from './home.module.css';
import ProjectBarCard from 'components/GeneMapper/projectBarCard';
import SearchIcon from '@mui/icons-material/Search';

const theme = createTheme({
  palette: {
    primary: {
      main: '#000000',
    },
    secondary: {
      main: '#5676E4',
    },
  },
  // typography: {
  //   fontFamily: 'Lato',
  //   fontSize: 16,
  //   // the size of save
  // },
});
// theme.typography.h3 = {
//   fontSize: '1.2rem',
//   '@media (min-width:600px)': {
//     fontSize: '1.5rem',
//   },
//   [theme.breakpoints.up('md')]: {
//     fontSize: '2rem',
//   },
// };
const themeIcon = createTheme({
  palette: {
    primary: {
      main: '#5676E4',
    },
  },
});
// dummy projects
function GeneMapperHome() {
  const [projects, setProjects] = useState([{
    id: 1,
    name: 'Leanne Graham',
    username: 'Bret',
    status: 'DONE',
  },
  {
    id: 2,
    name: 'Ervin Howell',
    username: 'Antonette',
    status: 'IN PROGRESS',
  },
  {
    id: 3,
    name: 'Max Musterman',
    username: 'Antonette',
    status: 'UPLOAD FAILED',
  },

  ]);
  const [findString, setFindString] = useState('');
  const [foundProjects, setFoundProjects] = useState({}, {});

  const filterFunktion = (e) => {
    const keyword = e.target.value;

    if (keyword !== '') {
      const results = projects.filter((project) => project.name.toLowerCase()
        .startsWith(keyword.toLowerCase()));
      setFoundProjects(results);
    } else {
      setFoundProjects(projects);
      // If the text field is empty, show all users
    }

    setFindString(keyword);
  };

  return (
    <div>
      <ThemeProvider theme={theme}>
        <Box
          sx={{
            marginLeft: '5%',
            width: '85%',
            display: 'flex',
            flexDirection: 'row',
            justifyContent: 'space-between',
            paddingBottom: '2em',
          }}
          className={styles.stack}
        >
          <Stack direction="row" className="stack">
            <Typography sx={{ fontWeight: 600, fontSize: '20px', marginTop: 0.5 }} className={styles.title}>Your Mappings </Typography>
            <ThemeProvider theme={themeIcon}>
              <PlusIcon className={styles.plus} />
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
            onChange={filterFunktion}
          />
        </Box>
        <div>
          {foundProjects && foundProjects.length > 0 ? (
            foundProjects.map((project) => (
              <ProjectBarCard name={project.name} status={project.status} />
            ))
          ) : (
            projects.map((project) => (
              <ProjectBarCard name={project.name} status={project.status} />
            ))
          )}
        </div>

      </ThemeProvider>
    </div>
  );
}

export default GeneMapperHome;
