/* eslint-disable */

import React from 'react'
import { Switch, Route, Redirect } from 'react-router-dom'
import { Box } from '@mui/material'
import Grid from '@mui/material/Grid';
import { useState, useEffect } from 'react'

import { TabGroup } from 'components/Tab'
import Search from 'components/Search'
import Filter from 'components/Filter'
import NavBar from 'components/NavBar'
import Breadcrumb from "components/Breadcrumb"
import AtlasCard from 'components/Cards/AtlasCard'
import DatasetCard from 'components/Cards/DatasetCard'
import { ModelCard } from 'components/Cards/ModelCard'
import Mapper from "components/Mapper"

import './Explore.css'

const tmpObj = [
    {
        label: "ATLASES",
        path: "/explore/atlases"
    },
    {
        label: "MODELS",
        path: "/explore/models"
    },
    {
        label: "DATASETS",
        path: "/explore/datasets"
    }
]

const atlasesGrid = (
    <Box className='atlasContainer'>
        <Grid container spacing={3}>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <AtlasCard width="300px" height="500px" title='Human - PBMC' imgLink='https://3-h.de/wp-content/uploads/grey-gradient-background.jpeg' modalities=' RNA, ADT' cellsInReference='161,764' species='Human' />
            </Grid>
        </Grid>
    </Box >
)

const modelsGrid = (
    <Box className='cardsContainer'>
        <Grid container spacing={3}>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <ModelCard title="Model 1" description="This is a short description" />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>

            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>

            </Grid>
        </Grid>
    </Box>
)

const datasetsGrid = (
    <Box className='cardsContainer'>
        <Grid container spacing={3}>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <DatasetCard title="title1" category="category1" />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <DatasetCard title="title2" category="category1" />
            </Grid>
            <Grid item xs={12} sm={6} md={4} lg={3}>
                <DatasetCard title="title3" category="category1" />
            </Grid>
        </Grid>
    </Box>
)

const Explore = () => {

    const [value, setValue] = useState(0)

    return (
        <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }} >

            <Box>
                <NavBar />
                <h1>NavBar goes here</h1>
            </Box>

            <Box sx={{ alignSelf: 'center', width: '65%', marginTop: '2%' }}>
                <Breadcrumb fontSize={1} />
            </Box>

            <Box sx={{ alignSelf: 'center', width: '65%', marginBlock: '2%' }}>
                <Search filterComponent={<Filter references={["test", "test"]} categories={["category1", "category2"]} />} handleSearch={(e) => console.log(e)} />
            </Box>

            <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'center', alignSelf: 'center', width: '80%' }}>
                {/* /explore/atlases */}
                <TabGroup value={value} setValue={setValue} tabsInfo={tmpObj} />
                <Switch>
                    <Route path="/explore/atlases" render={() => atlasesGrid} />
                    <Route path="/explore/models" render={() => modelsGrid} />
                    <Route path="/explore/datasets" render={() => datasetsGrid} />
                    <Redirect to="/explore/atlases" />
                </Switch>

            </Box>

            <Mapper />
        </Box>
    )
}

export default Explore