import React, { useCallback, useState } from 'react';
import TextField from '@mui/material/TextField';
import Stack from '@mui/material/Stack';
import Button from '@mui/material/Button';
import Select from '@mui/material/Select';
import FormControl from '@mui/material/FormControl';
import MenuItem from '@mui/material/MenuItem';
import InputLabel from '@mui/material/InputLabel';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import * as IoIcons from 'react-icons/io';
import * as CGIcons from 'react-icons/cg';
import completeImg from '../../../assets/complete.png';
import styles from './popup.module.css';
import projects from '../Sidebar/SidebarData';

const theme = createTheme({
  palette: {
    primary: {
      main: '#1888ff',
    },
    secondary: {
      main: '#6B7379',
    },
  },
});

function PopUp() {
  // Booleans
  const [showPopup, setShowPopup] = React.useState(true);

  if (showPopup) {
    return (
      <>
        <div className={styles.background} />
        <PopUpContent setShowPopup={setShowPopup} />
      </>
    );
  }
  return null;
}

function PopUpContent(props) {
  // make projects stateful
  const [projectList, setProjects] = useState(projects);

  // Booleans
  // eslint-disable-next-line no-unused-vars
  const [newProject, setNewProject] = useState(false);

  // Annotation
  const [annotationDetails, setAnnotationDetails] = useState({
    annotationName: '',
    newProjectName: '',
    projectName: '',
  });

  const handleTextChange = useCallback((e) => {
    setAnnotationDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [annotationDetails]);

  function addProject(pName) {
    // create new project
    console.log('PNAME: ');
    const newProj = {
      name: pName,
      path: `/${pName}`,
      iconClosed: <IoIcons.IoIosArrowForward />,
      iconOpened: <IoIcons.IoIosArrowDown />,
      subNav: [],
    };
    // add new project
    setProjects([...projectList, newProj]);
  }

  function addAnnotation(pName) {
    const newAnnot = {
      name: annotationDetails.annotationName,
      path: `/overview/${annotationDetails.annotationName}`,
      icon: <CGIcons.CgFileDocument />,
    };
    const idx = projectList.findIndex((e) => e.name === pName);
    if (idx >= 0) {
      setProjects(projectList[idx].subNav.concat({ newAnnot }));
    }
  }

  // Button: Add Annotation
  const addNewCelltypeAnnotation = useCallback(() => {
    const projName = newProject && annotationDetails.newProjectName !== ''
      ? annotationDetails.newProjectName : annotationDetails.projectName;

    if (projectList.filter((e) => e.name === projName).length === 0) {
      addProject(projName);
    }
    addAnnotation(projName);
    props.setShowPopup(false);
  }, [projectList]);

  return (
    <div className={styles.popUpBase}>

      <div className={styles.popUpHeadline}>
        Processing complete!
        <img
          alt="complete"
          src={completeImg}
          style={{ height: '40px', paddingLeft: '20px' }}
        />
      </div>

      <div style={{ paddingBlock: '40px', paddingInline: '100px' }}>
        <div className={styles.divider} />
      </div>

      <Stack
        spacing={5}
        className={styles.inputComponentList}
        style={{ left: '50px' }}
      >

        <div className={styles.inputComponent}>
          <div className={styles.popUpInputText}>
            New Cell Type Annotation
          </div>

          <TextField
            required
            id="annotationName"
            label="Annotation Name"
            type="text"
            style={{ width: '350px' }}
            onChange={handleTextChange}
          />
        </div>

        <div className={styles.inputComponent}>
          <div className={styles.popUpInputText}>
            Associated Project
          </div>

          <TextField
            required
            id="newProjectName"
            label="New Project Name"
            type="text"
            style={{ width: '350px' }}
            onChange={handleTextChange}
          />
        </div>

        <div>
          <FormControl sx={{ width: 350 }}>
            <InputLabel
              id="demo-simple-select-autowidth-label"
              style={{ textAlign: 'left' }}
            >
              Project Name*
            </InputLabel>
            <Select
              labelId="demo-simple-select-autowidth-label"
              id="projectName"
              autoWidth
              label="Project Name*"
              onChange={handleTextChange}
            >
              {/* Elements to choose from */}
              <div style={{ width: '350px' }} />
              { projectList.map((item, index) => (
                <MenuItem value={index}>
                  {item.name}
                </MenuItem>
              ))}

            </Select>
          </FormControl>
        </div>

        <ThemeProvider theme={theme}>
          <Stack
            spacing={2}
            direction="row"
            style={{ height: '55px', marginTop: '50px' }}
          >
            <Button
              variant="outlined"
              color="secondary"
              style={{ width: '150px', borderRadius: '10px' }}
            >
              Cancel
            </Button>

            <Button
              variant="contained"
              color="primary"
              style={{ width: '185px', borderRadius: '10px', fontWeight: 'bold' }}
              disabled={annotationDetails.annotationName === ''
                     || ((annotationDetails.newProjectName === '' && !newProject)
                     || (annotationDetails.projectName === '' && newProject))}
              onClick={addNewCelltypeAnnotation}
            >
              Add Annotation
            </Button>
          </Stack>
        </ThemeProvider>
      </Stack>
    </div>
  );
}

export default PopUp;
