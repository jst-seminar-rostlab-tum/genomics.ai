import React, { useState, useEffect, useCallback } from 'react';
import {
  Switch,
  Route,
  Redirect,
  useHistory,
  useLocation,
} from 'react-router-dom';
import { Box } from '@mui/material';
import Grid from '@mui/material/Grid';
import { TabGroup } from 'components/Tab';
import Search from 'components/Search';
import { Filter } from 'components/Filter/Filter';
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
import LearnMore from './LearnMore';

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
  const [searchValue, setSearchValue] = useState('');

  const handleSearch = (e) => {
    setSearchValue(e);
    console.log(e);
  };

  const [atlases, setAtlases] = useState([]);
  const [models, setModels] = useState([]);

  const path = useLocation().pathname;

  useEffect(() => {
    AtlasService.getAtlases()
      .then((atlases) => setAtlases(atlases))
      .catch((err) => console.log(err));

    ModelsService.getModels()
      .then((models) => setModels(models))
      .catch((err) => console.log(err));
  }, []);

  useEffect(() => {
    if (selectedAtlas || selectedModel) setMapperVisible(true);
    if (!selectedAtlas && !selectedModel) setMapperVisible(false);
  }, [selectedAtlas, selectedModel]);

  const AtlasesGrid = ({ atlases, path }) => (
    <Box className="atlasContainer" sx={{ height: '70vh' }}>
      <Grid container spacing={3}>
        {atlases && atlases.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => setSelectedAtlas(atlas)}
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/${atlas._id}`}
            />
          </Grid>
        ))}
        {atlases && atlases.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => setSelectedAtlas(atlas)}
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/${atlas._id}`}
            />
          </Grid>
        ))}
        {atlases && atlases.map((atlas) => (
          <Grid key={atlas._id} item xs={12} sm={6} md={4} lg={3}>
            <AtlasCard
              onClick={() => setSelectedAtlas(atlas)}
              atlasId={atlas._id}
              imgLink={atlas.previewPictureURL}
              species={atlas.species}
              modalities={atlas.modalities}
              title={atlas.name}
              learnMoreLink={`${path}/${atlas._id}`}
            />
          </Grid>
        ))}
      </Grid>
    </Box>
  );

  const ModelsGrid = ({ models, path }) => (
    <Box className="cardsContainer">
      <Grid container spacing={3}>
        {models && models.map((model) => (
          <Grid key={model._id} item xs={12} sm={6} md={4} lg={3}>
            <ModelCard
              onClick={() => setSelectedModel(model)}
              title={`Model ${model.name}`}
              description={model.description}
              learnMoreLink={`${path}/${model._id}`}
            />
          </Grid>
        ))}
      </Grid>
    </Box>
  );
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
  const tabMenu = (props) => (
    <>
      <Box sx={{ alignSelf: 'center', width: '60%', marginBlock: '2%' }}>
        <Search
          filterComponent={(
            <Filter
              references={['test', 'test']}
              categories={['category1', 'category2']}
            />
)}
          handleSearch={handleSearch}
          value={searchValue}
        />
      </Box>
      <TabGroup value={value} setValue={setValue} tabsInfo={tmpObj} />
      {value === 0 ? <AtlasesGrid atlases={props.atlases} path={props.path} /> : null }
      {value === 1 ? <ModelsGrid models={props.models} path={props.path} /> : null }
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
        />
      </Box>

      <Box sx={{ alignSelf: 'center', width: '60%', marginTop: '2%' }}>
        <Breadcrumb fontSize={1} actions={{ explore: () => setValue(0) }} />
      </Box>

      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          alignSelf: 'center',
          width: '60%',
        }}
      >
        {/* /explore/atlases */}
        <Switch>
          <Route
            exact
            path="/explore/atlases"
            render={() => atlases && tabMenu({ atlases, path })}
          />
          <Route
            path="/explore/models"
            render={() => models && tabMenu({ models, path })}
          />
          <Route path="/explore/datasets" render={() => atlases && tabMenu(models, path, 'datasets')} />
          <Route path="/explore/atlases/:id" render={() => <LearnMore setSelectedAtlas={setSelectedAtlas} setSelectedModel={setSelectedModel} />} />
          <Redirect to="/explore/atlases" />
        </Switch>
      </Box>

      <Mapper
        mapperAtlas={selectedAtlas ? selectedAtlas.name : null}
        mapperModel={selectedModel ? selectedModel.name : null}
        setSelectedAtlas={setSelectedAtlas}
        setSelectedModel={setSelectedModel}
        open={mapperVisible}
        fabOnClick={() => setMapperVisible(!mapperVisible)}
      />
    </Box>
  );
};

export default Explore;
