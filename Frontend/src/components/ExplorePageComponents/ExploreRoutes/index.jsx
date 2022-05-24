import React from 'react';
import {
  Switch, Route, Redirect,
} from 'react-router-dom';

import AtlasResult from 'views/Explore/AtlasResult';
import LearnMoreAtlas from 'views/Explore/LearnMoreAtlas';
import LearnMoreModel from 'views/Explore/LearnMoreModel';
import GeneMapperState from 'views/GeneMapper/GeneMapperState';
import GeneMapperHome from 'views/GeneMapper/Home';

const ExploreRoutes = ({
  path, atlases, models, selectedAtlas, selectedModel,
}) => (
  <Switch>
    <Route
      exact
      path={`${path}/atlases`}
      render={() => atlases}
    />
    <Route
      exact
      path={`${path}/models`}
      render={() => models}
    />
    <Route exact path={`${path}/models/:id`} render={() => <LearnMoreModel />} />
    <Route exact path={`${path}/atlases/:id/visualization`} render={() => <AtlasResult />} />
    <Route exact path={`${path}/atlases/:id`} render={() => <LearnMoreAtlas />} />
    <Route exact path={`${path}/genemapper/create`}>
      <GeneMapperState
        path={path}
        defaultSelectedAtlas={selectedAtlas}
        defaultSelectedModel={selectedModel}
      />
    </Route>
    <Route exact path={`${path}/genemapper`}>
      <GeneMapperHome />
    </Route>
    <Redirect to={`${path}/atlases`} />
  </Switch>
);

export default ExploreRoutes;
