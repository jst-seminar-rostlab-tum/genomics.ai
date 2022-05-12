import React, { useState, useEffect, useCallback } from 'react';
import {
  Switch, Route, Redirect, useHistory, useLocation, useParams,
  useRouteMatch, Link
} from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { TabGroup } from 'components/Tab';
import Search from 'components/Search';
import Filter from 'components/ExplorePageComponents/Filter';
import NavBar from 'components/NavBar';
import Breadcrumb from 'components/Breadcrumb';
import AtlasCard from 'components/Cards/AtlasCard';
import { ModelCard } from 'components/Cards/ModelCard';
import Mapper from 'components/Mapper';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import './Explore.css';

import DatasetCard from 'components/Cards/DatasetCard';
import ModelsService from 'shared/services/Models.service';
import AtlasService from 'shared/services/Atlas.service';
import LearnMoreAtlas from './LearnMoreAtlas';
import LearnMoreModel from './LearnMoreModel';
import ResultVisualization from 'components/GeneMapper/ResultVisualization';
import AtlasResult from './AtlasResult';

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
    label: 'DATASETS',
    path: '/explore/datasets',
  },
];

const datasets = [
  {
    title: 'dataset 1',
    category: 'human',
  },
  {
    title: 'dataset 2',
    category: 'animal',
  },
];

const Explore = () => {
  const [value, setValue] = useState(0);
  const [selectedAtlas, setSelectedAtlas] = useState(null);
  const [selectedModel, setSelectedModel] = useState(null);
  const [mapperVisible, setMapperVisible] = useState(false);
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);
  const { search } = useLocation();
  const searchParams = new URLSearchParams(search);
  const searchedKeyword = searchParams.get('keyword') || '';
  const { path } = useRouteMatch();
  const history = useHistory();
  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);

  // function to update the state in the URL
  const updateQueryParams = (param, value) => {
    const params = new URLSearchParams(history.location.search);
    if (value) {
      params.set(param, value);
    } else {
      params.delete(param);
    }
    history.push({
      pathname: history.location.pathname,
      search: params.toString(),
    });
  };

  useEffect(() => {
    AtlasService.getAtlases()
      .then((newAtlases) => setAtlases(newAtlases))
      .catch((err) => console.log(err));

    ModelsService.getModels()
      .then((newModels) => setModels(newModels))
      .catch((err) => console.log(err));
  }, []);

  const searchedKeywordChangeHandler = (value) => {
    updateQueryParams('keyword', value);
  };

  useEffect(() => {
    if (selectedAtlas || selectedModel) setMapperVisible(true);
    if (!selectedAtlas && !selectedModel) setMapperVisible(false);
  }, [selectedAtlas, selectedModel]);

  function applyAtlasFilters() {
    const searchedAtlases = atlases.filter(
      (item) => item.name.toLowerCase().includes(searchedKeyword.toLowerCase()),
    );
    if (searchParams.get('sortBy') === 'name' || searchParams.get('sortBy') === null) {
      searchedAtlases.sort((a, b) => {
        const nameA = a.name.toUpperCase();
        const nameB = b.name.toUpperCase();
        if (nameA < nameB) {
          return -1;
        }
        if (nameA > nameB) {
          return 1;
        }
        return 0;
      });
    } else if (searchParams.get('sortBy') === 'numberOfCells') {
      searchedAtlases.sort((a, b) => a.numberOfCells - b.numberOfCells);
    }
    return searchedAtlases;
  }

  function applyModelFilters() {
    const searchedModels = models.filter(
      (item) => item.name.toLowerCase().includes(searchedKeyword.toLowerCase()),
    );
    if (searchParams.get('sortBy') === 'name' || searchParams.get('sortBy') === null) {
      searchedModels.sort((a, b) => {
        const nameA = a.name.toUpperCase();
        const nameB = b.name.toUpperCase();
        if (nameA < nameB) {
          return -1;
        }
        if (nameA > nameB) {
          return 1;
        }
        return 0;
      });
    }
    return searchedModels;
  }

  const AtlasesGrid = () => {
    const atlasesFiltered = applyAtlasFilters();
    return (
      <Box className="atlasContainer" sx={{ height: '70vh' }}>
        <Grid container spacing={3}>
          {atlasesFiltered && atlasesFiltered.map((atlas) => (
            <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
              <AtlasCard
                onClick={() => history.push(`${path}/atlases/${atlas._id}/visualization`)}
                atlasId={atlas._id}
                imgLink={atlas.previewPictureURL}
                species={atlas.species}
                modalities={atlas.modalities}
                title={atlas.name}
                learnMoreLink={`${path}/atlases/${atlas._id}`}
              />
            </Grid>
          ))}
          {atlasesFiltered && atlasesFiltered.map((atlas) => (
            <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
              <AtlasCard
                onClick={() => setSelectedAtlas(atlas)}
                atlasId={atlas._id}
                imgLink={atlas.previewPictureURL}
                species={atlas.species}
                modalities={atlas.modalities}
                title={atlas.name}
                learnMoreLink={`${path}/atlases/${atlas._id}`}
              />
            </Grid>
          ))}
          {atlasesFiltered && atlasesFiltered.map((atlas) => (
            <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
              <AtlasCard
                onClick={() => setSelectedAtlas(atlas)}
                atlasId={atlas._id}
                imgLink={atlas.previewPictureURL}
                species={atlas.species}
                modalities={atlas.modalities}
                title={atlas.name}
                learnMoreLink={`${path}/atlases/${atlas._id}`}
              />
            </Grid>
          ))}
        </Grid>
      </Box>
    );
  };

  const ModelsGrid = () => {
    const modelsFiltered = applyModelFilters();
    return (
      <Box className="cardsContainer">
        <Grid container spacing={3}>
          {modelsFiltered && modelsFiltered.map((model) => (
            <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
              <ModelCard
                onClick={() => setSelectedModel(model)}
                title={`Model ${model.name}`}
                description={model.description}
                learnMoreLink={`${path}/models/${model._id}`}
              />
            </Grid>
          ))}
        </Grid>
      </Box>
    );
  };
  const DataSetGrids = () => (
    <Box className="cardsContainer">
      <Grid container spacing={3}>
        {datasets.map((d, index) => (
          <Grid key={index} item xs={12} sm={6} md={4} lg={3}>
            <DatasetCard title={d.title} category={d.category} />
          </Grid>
        ))}
      </Grid>
    </Box>
  );
  const tabMenu = () => (
    <>
      <TabGroup value={value} setValue={setValue} tabsInfo={tmpObj} />
      {value === 0 ? <AtlasesGrid /> : null }
      {value === 1 ? <ModelsGrid /> : null }
      {value === 2 ? <DataSetGrids /> : null }
    </>

  );

  const onLoginClicked = useCallback(() => {
    console.log('login');
    setRegistrationFormVisible(false);
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

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

  const tmp_elems = useLocation().pathname.split('/');
  const elems = tmp_elems.map((elem, index) => {
    if (index === 3) {
      if (tmp_elems[2] === 'atlases') return atlases.filter((x) => x._id === elem)[0] ? atlases.filter((x) => x._id === elem)[0].name : elem;
      if (tmp_elems[2] === 'models') return models.filter((x) => x._id === elem)[0] ? models.filter((x) => x._id === elem)[0].name : elem;
    }
    return elem;
  });

  const executeScroll = () => history.push({pathname: '/', state: {contact_us: true}})

  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center',
      }}
    >
      {isLoginFormVisible && (
        <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />
      )}
      {isRegistrationFormVisible && (
        <RegistrationForm
          visible={isRegistrationFormVisible}
          onClose={onRegistrationFormClosed}
        />
      )}

      <Box>
        <NavBar
          position="relative"
          onLoginClicked={onLoginClicked}
          onSignUpClicked={onSignUpClicked}
          executeScroll={executeScroll}
        />
      </Box>

      <Box sx={{ alignSelf: 'center', width: '60%', marginTop: '2%' }}>
        <Breadcrumb elems={elems} fontSize={1} actions={{ explore: () => setValue(0) }} />
      </Box>

      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          alignSelf: 'center',
          width: '60%',
        }}
      >
        <Box sx={{ alignSelf: 'center', width: '100%', marginBlock: '2%' }}>
          <Search
            filterComponent={(
              <Filter
                searchParams={searchParams}
                updateQueryParams={updateQueryParams}
                path={path}
              />
            )}
            handleSearch={searchedKeywordChangeHandler}
            value={searchedKeyword}
          />
        </Box>
        {/* /explore/atlases */}
        <Switch>
          <Route
            exact
            path="/explore/atlases"
            render={() => atlases && tabMenu()}
          />
          <Route
            exact
            path="/explore/models"
            render={() => models && tabMenu()}
          />
          <Route exact path="/explore/models/:id" render={() => <LearnMoreModel />} />
          <Route exact path="/explore/datasets" render={() => atlases && tabMenu()} />
          <Route exact path="/explore/atlases/:id/visualization" render={() => <AtlasResult />} />
          <Route exact path="/explore/atlases/:id" render={() => <LearnMoreAtlas />} />
          <Redirect to="/explore/atlases" />
        </Switch>
      </Box>

      {/* NOT NEEDED FOR NOW */}
      {/* <Mapper
        mapperAtlas={selectedAtlas ? selectedAtlas.name : null}
        mapperModel={selectedModel ? selectedModel.name : null}
        setSelectedAtlas={setSelectedAtlas}
        setSelectedModel={setSelectedModel}
        open={mapperVisible}
        fabOnClick={() => setMapperVisible(!mapperVisible)}
      /> */}
    </Box>
  );
};

export default Explore;
