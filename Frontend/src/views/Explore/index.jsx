import React, { useState, useEffect, useCallback } from 'react';
import {
  Switch, Route, Redirect, useHistory,
} from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { TabGroup } from 'components/Tab';
import Search from 'components/Search';
import { Filter } from 'components/Filter/Filter';
import NavBar from 'components/NavBar';
import Breadcrumb from 'components/Breadcrumb';
import AtlasCard from 'components/Cards/AtlasCard';
import DatasetCard from 'components/Cards/DatasetCard';
import { ModelCard } from 'components/Cards/ModelCard';
import Mapper from 'components/Mapper';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import './Explore.css';
import { useState, useEffect, useCallback } from 'react'

import ModelsService from 'shared/services/Models.service'
import AtlasService from 'shared/services/Atlas.service'

import { useLocation } from 'react-router-dom';
import { useAuth } from 'shared/context/authContext';

const tmpObj = [
  {
    label: 'ATLASES',
    path: '/explore/atlases',
  },
  {
    label: 'MODELS',
    path: '/explore/models',
  },
  {
    label: "DATASETS",
    path: "/explore/datasets"
  }
]

const atlasesGrid = (atlases, path) => (
  <Box className='atlasContainer'>
    <Grid container spacing={3}>
      {
        atlases.map((atlas) => (
          <Grid item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard width="300px" height="500px" imgLink={atlas.previewPictureURL} species={atlas.species} modalities={atlas.modalities} title={atlas.title} learnMoreLink={`${path}/${atlases._id}`} />
          </Grid>
        ))
      }
      <Grid item xs={12} sm={6} md={4} lg={3}>
        <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
      </Grid>
    </Grid>
  </Box >
)

const modelsGrid = (models, path) => (
  <Box className='cardsContainer'>
    <Grid container spacing={3}>
      {
        models.map((model) => (
          <Grid item xs={12} sm={6} md={4} lg={3}>
            <ModelCard title={`Model ${model.name}`} description={model.description} learnMoreLink={`${path}/${model._id}`} />
          </Grid>
        ))
      }
    </Grid>
  </Box>
)

const datasetsGrid = (
  <Box className="cardsContainer">
    <Grid container spacing={3}>
      {datasetTitles.map((title, index) => (
        <Grid item xs={12} sm={6} md={4} lg={3}>
          <DatasetCard title={title} category={datasetCategories[index]} />
        </Grid>
      ))}
    </Grid>
  </Box>
)

const Explore = () => {

  const [, setUser] = useAuth()

  const Explore = () => {
    const [value, setValue] = useState(0);
    const [selectedAtlas, setSelectedAtlas] = useState(-1);
    const [selectedModel, setSelectedModel] = useState(-1);
    const [mapperVisible, setMapperVisible] = useState(false);
    const [isLoginFormVisible, setLoginFormVisible] = useState(false);
    const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);
    const [searchValue, setSearchValue] = useState('');
    const history = useHistory();

    const handleSearch = (e) => {
      setSearchValue(e);
      console.log(e);
    };

    const handleAtlasMapClick = (index) => {
      setSelectedAtlas(index);
      if (!mapperVisible) {
        setMapperVisible(true);
        if (selectedModel === -1) {
          history.push('/explore/models');
        }
      }
    };

    const [atlases, setAtlases] = useState([])
    const [models, setModels] = useState([])

    const path = useLocation().pathname

    console.log(path)

    useEffect(() => {
      AtlasService.getAtlases()
        .then((atlases) => setAtlases(atlases))
        .catch((err) => console.log(err))

      ModelsService.getModels()
        .then((models) => setModels(models))
        .catch((err) => console.log(err))

      console.log(atlases)
      console.log(models)
    }, [])

    const onLoginClicked = useCallback(() => {
      console.log("login")
      setRegistrationFormVisible(false)
      setLoginFormVisible(true)
    }, [setLoginFormVisible])

    const onSignUpClicked = useCallback(() => {
      console.log('register');
      setLoginFormVisible(false);
      setRegistrationFormVisible(true);
    }, [setRegistrationFormVisible]);

    const onLoginFormClosed = useCallback(() => {
      setLoginFormVisible(false);
    }, [setLoginFormVisible]);

    const onRegistrationFormClosed = useCallback(() => {
      setRegistrationFormVisible(false);
    }, [setRegistrationFormVisible]);

    const handleModelMapClick = (index) => {
      setSelectedModel(index);
      if (!mapperVisible) {
        setMapperVisible(true);
      }
      if (selectedAtlas === -1) {
        history.push('/explore/atlases');
      }
    };

    const atlasesGrid = (
      <Box className="atlasContainer">
        <Grid container spacing={3}>
          {atlasTitles.map((title, index) => (
            <Grid item xs={12} sm={6} md={4} lg={3}>
              <AtlasCard width="270px" height="500px" index={index} title={title} imgLink={atlasImgLinks[index]} modalities={atlasModalities[index]} cellsInReference={atlasCellsInReference[index]} species={atlasSpecies[index]} onClick={handleAtlasMapClick} />
            </Grid>
          ))}
        </Grid>
      </Box>
    );

    const modelsGrid = (
      <Box className="cardsContainer">
        <Grid container spacing={3}>
          {modelTitles.map((title, index) => (
            <Grid item xs={12} sm={6} md={4} lg={3}>
              <ModelCard
                title={title}
                description={modelDescriptions[index]}
                onClick={handleModelMapClick}
                index={index}
              />
            </Grid>
          ))}
        </Grid>
      </Box>
    );

    return (
      <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>

        {isLoginFormVisible
          && <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />}
        {isRegistrationFormVisible
          && (
            <RegistrationForm
              visible={isRegistrationFormVisible}
              onClose={onRegistrationFormClosed}
            />
          )}

        <Box>
          <NavBar position="relative" onLoginClicked={onLoginClicked} onSignUpClicked={onSignUpClicked} />
        </Box>

        <Box sx={{ alignSelf: 'center', width: '65%', marginTop: '2%' }}>
          <Breadcrumb fontSize={1} actions={{ explore: () => setValue(0) }} />
        </Box>

        <Box sx={{ alignSelf: 'center', width: '65%', marginBlock: '2%' }}>
          <Search filterComponent={<Filter references={['test', 'test']} categories={['category1', 'category2']} />} handleSearch={handleSearch} value={searchValue} />
        </Box>

        <Box sx={{
          display: 'flex', flexDirection: 'column', justifyContent: 'center', alignSelf: 'center', width: '80%',
        }}
        >
          {/* /explore/atlases */}
          <TabGroup value={value} setValue={setValue} tabsInfo={tmpObj} />
          <Switch>
            <Route path="/explore/atlases" render={() => atlasesGrid(atlases, path)} />
            <Route path="/explore/models" render={() => modelsGrid(models, path)} />
            <Route path="/explore/datasets" render={() => datasetsGrid} />
            <Redirect to="/explore/atlases" />
          </Switch>

        </Box>

        <Mapper
          mapperAtlas={atlasTitles[selectedAtlas] ?? null}
          mapperModel={modelTitles[selectedModel] ?? null}
          setSelectedAtlas={setSelectedAtlas}
          setSelectedModel={setSelectedModel}
          open={mapperVisible}
          fabOnClick={() => setMapperVisible(!mapperVisible)}
        />
      </Box>
    );
  };
}

export default Explore;
