"use strict";angular.module("inspinia",["ngAnimate","ngCookies","ngTouch","ngSanitize","ngResource","ui.router","ui.bootstrap","highcharts-ng","ui.select","angular-venn","angularMoment"]).config(["$stateProvider","$urlRouterProvider",function(t,e){t.state("index",{"abstract":!0,url:"/index",templateUrl:"components/common/content.html"}).state("index.main",{url:"/main",templateUrl:"app/main/main.html",data:{pageTitle:"Dashboard"}}).state("index.school",{url:"/school",templateUrl:"components/school/school.html",data:{pageTitle:"School Analytics"}}),e.otherwise("/index/main")}]),angular.module("inspinia").controller("TabsDemoCtrl",["$scope","$rootScope","$window",function(t){t.tabs=[{title:"Global",content:"Dynamic content 1"},{title:"China",content:"Dynamic content 2",disabled:!0}],t.click=function(){setTimeout(function(){window.dispatchEvent(new Event("resize"))},100)}}]),angular.module("inspinia").controller("UserRankingCrol",["$scope","$rootScope","$http",function(t,e,n){var i=[{name:"FX",id:"1"},{name:"FC",id:"2"},{name:"SC",id:"3"}];t.rankListLoaded=!1;var a=[];t.$watch(function(){return e.selectedSchool},function(){e.selectedSchool&&(t.products="CN"==e.regionSelected?i:[i[0]],t.selected=t.products[0],t.loadSchoolData(e.selectedSchool))},!0),t.loadSchoolData=function(i){t.loaded=!1,t.userRankList=[],t.rankListLoaded=!1,n.get("getUserRanking?school_key="+i+"&app="+t.selected.name+"&db="+e.regionSelected).success(function(e){a=e,t.rankListLoaded=!0,t.userRankList=e.slice(0,10)})},t.loadMore=function(){t.userRankList=a},t.changeProduct=function(){t.loadSchoolData(e.selectedSchool)}}]),angular.module("inspinia").controller("SchoolCtrl",["$scope","$rootScope","$http",function(t,e,n){t.regions=[],t.schoolsList=[],t.schools=[],t.school={},t.compareList=[],n.get("getSchools").success(function(e){for(var n in e)e[n].school_name&&(e[n].school_region&&t.regions.indexOf(e[n].school_region)<0&&t.regions.push(e[n].school_region),"INOTHER"!=e[n].school_key&&t.schoolsList.push(e[n]))}),t.addSchool=function(){t.school.selected="Please select a school to compare ...",t.isCompard=!0},t.loadSchoolData=function(e){n.get("getSchoolUserGrowth?school_key="+e.school_key+"&db="+e.school_region).success(function(n){t.isCompard||(t.chartConfig.series=[]);for(var i=[],a=0;a<n.length;a++)i.push([Date.parse(n[a].date),n[a].total_user_num]);t.chartConfig.series.push({type:"spline",name:e.school_name,data:i})})},t.onSearch=function(){},t.loadCountData=function(e){t.totalUser=0,t.monthUser=0,t.weekUser=0,n.get("getSchoolUserCount?school_key="+e.school_key+"&db="+e.school_region).success(function(e){var n=e;t.totalUser=n.totalUser,t.monthUser=n.monthlyActiveUser,t.weekUser=n.weeklyActiveUser})},t.onSchoolSelect=function(n){t.count_show=!0,t.loadSchoolData(n),t.loadCountData(n),e.selectedSchool=n.school_key},t.onRegionSelect=function(n){e.regionSelected=n,t.school_show=!0,t.school.selected=void 0,t.schools=[];for(var i in t.schoolsList)t.schoolsList[i].school_region&&t.schoolsList[i].school_region==n&&t.schools.push(t.schoolsList[i])},t.chartConfig={chart:{backgroundColor:null,zoomType:"x"},title:{text:"User Growth"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("ContestCrol",["$scope","$rootScope","$http",function(t,e,n){var i=[{name:"FX",id:"1"},{name:"FC",id:"2"},{name:"SC",id:"3"}];t.$watch(function(){return e.selectedSchool},function(){e.selectedSchool&&(t.products="CN"==e.regionSelected?i:[i[0]],t.selected=t.products[0],t.loadContentData(e.selectedSchool))},!0),t.loadContentData=function(i){t.contestList=[],n.get("getContestByApp?school_key="+i+"&app="+t.selected.name+"&db="+e.regionSelected).success(function(e){t.contestList=e})},t.changeProduct=function(){t.loadContentData(e.selectedSchool)}}]);var app=angular.module("inspinia");app.controller("MapCtrl",["$scope","$location","$http",function(t,e,n){t.data=function(){n.get("getRegionDist?file=Global_Country_TU").success(function(e){e=e.data.Global_Country_TU;var n={};t.markers=[];for(var i=0;i<e.length;i++)e[i].Country&&(n[e[i].Country]=parseInt(e[i]["Total Users"])),"HK"==e[i].Country&&t.markers.push({latLng:[22.15,114.1],name:"Hong Kong",values:parseInt(e[i]["Total Users"])});t.datamap=n})}}]),app.directive("map",function(){return{restrict:"EAC",link:function(t,e){var n=null;t.$watch("datamap",function(){console.log(n),t.datamap&&($(e).width("auto"),$(e).height(400),n=$(e).vectorMap({onRegionClick:function(t,e){console.info(e)},map:"world_mill_en",backgroundColor:"transparent",regionStyle:{initial:{fill:"#e4e4e4","fill-opacity":1,stroke:"none","stroke-width":0,"stroke-opacity":0}},series:{regions:[{values:t.datamap,scale:["#1ab394","#ED5565"],normalizeFunction:"polynomial"}]},onMarkerTipShow:function(e,n,i){n.html(n.html()+": "+(t.markers[i]?t.markers[i].values:""))},onRegionTipShow:function(e,n,i){n.html(n.html()+": "+(t.datamap[i]?t.datamap[i]:""))},markers:t.markers,markerStyle:{initial:{fill:"#F8E23B",stroke:"#383f47"}}}))})}}}),angular.module("inspinia").controller("VeenCtrl",["$scope","$http","$rootScope",function(t,e,n){t.loadData=function(i){e.get("getDataFromFile?file="+i).success(function(e){e=e.data;var i=[{sets:["China_FX"],size:e.China_FX_FC_SC_TU[0].China_FX},{sets:["China_FC"],size:e.China_FX_FC_SC_TU[0].China_FC},{sets:["China_SC"],size:e.China_FX_FC_SC_TU[0].China_SC},{sets:["China_FX","China_FC"],size:e.China_FX_FC_SC_TU[0].China_FX_FC},{sets:["China_FX","China_SC"],size:e.China_FX_FC_SC_TU[0].China_FX_SC},{sets:["China_FC","China_SC"],size:e.China_FX_FC_SC_TU[0].China_FC_SC},{sets:["China_FX","China_FC","China_SC"],size:e.China_FX_FC_SC_TU[0].China_FX_FC_SC}];n.overlapUsers=e.China_FX_FC_SC_TU[0].China_FX_FC_SC,t.vennData=i,setTimeout(function(){var t=d3.select("#venn"),e=d3.select("body").append("div").attr("class","venntooltip");t.selectAll("g").on("mouseover",function(n){venn.sortAreas(t,n),e.transition().duration(400).style("opacity",.9),e.text(n.size+" users");var i=d3.select(this).transition("tooltip").duration(400);i.select("path").style("stroke-width",3).style("fill-opacity",1==n.sets.length?.4:.1).style("stroke-opacity",1)}).on("mousemove",function(){e.style("left",d3.event.pageX+"px").style("top",d3.event.pageY-28+"px")}).on("mouseout",function(t){e.transition().duration(400).style("opacity",0);var n=d3.select(this).transition("tooltip").duration(400);n.select("path").style("stroke-width",0).style("fill-opacity",1==t.sets.length?.25:0).style("stroke-opacity",0)})},300)})}}]),angular.module("inspinia").controller("ChartCtrl",["$scope","$http","$rootScope",function(t,e,n){t.$on("changeSelector",function(t,e){console.info(e)}),t.loadData=function(i){e.get("getDataFromFile?file="+i).success(function(e){e=e.data;var a=i.split("|");t.chartConfig.series=[];for(var o=0;o<a.length;o++){for(var s=[],l=0;l<e[a[o]].length;l++)s.push([Date.parse(e[a[o]][l].date),parseInt(e[a[o]][l][t.valueName])]);t.chartConfig.series.push({type:"area",name:t[a[o]],data:s})}if("Date_TU_Global|Date_TU_China"==i){var r=e.Date_TU_Global,c=e.Date_TU_China;n.totalUsersChina=parseInt(c[c.length-1][t.valueName]),n.totalUsersGlobal=parseInt(r[r.length-1][t.valueName])}})},t.chartConfig={colors:["#A3E1D4","#ff0066","#2b908f","#90ee7e","#eeaaee","#55BF3B","#DF5353","#7798BF","#aaeeee"],chart:{backgroundColor:null,zoomType:"x"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("ChartComLineCtrl",["$scope","$http",function(t,e){t.loadData=function(n){e.get("getDataFromFile?file="+n).success(function(e){e=e.data;var n=[],i=[],a={};a={Others:0};for(var o={},s=0;s<e.Global_Country_TU.length&&4>s;s++)n.push(e.Global_Country_TU[s]),i.push(e.Global_Country_TU[s].Country),a[e.Global_Country_TU[s].Country]=0;for(var s=0;s<e.Date_Global_Country_TU.length;s++)if(i.indexOf(e.Date_Global_Country_TU[s].country)>=0){"TW"==e.Date_Global_Country_TU[s].country&&console.info(e.Date_Global_Country_TU[s]);var l=a[e.Date_Global_Country_TU[s].country]+parseInt(e.Date_Global_Country_TU[s]["Total Users"]);a[e.Date_Global_Country_TU[s].country]=l;var r=[Date.parse(e.Date_Global_Country_TU[s].date),l];o[e.Date_Global_Country_TU[s].country]&&o[e.Date_Global_Country_TU[s].country].length>0?o[e.Date_Global_Country_TU[s].country].push(r):(o[e.Date_Global_Country_TU[s].country]=[],o[e.Date_Global_Country_TU[s].country].push(r))}else{var l=a.Others+parseInt(e.Date_Global_Country_TU[s]["Total Users"]);a.Others=l;var r=[Date.parse(e.Date_Global_Country_TU[s].date),l];o.Others&&o.Others.length>0?o.Others.push(r):(o.Others=[],o.Others.push(r))}var c=[];i.push("Others");for(var h in i)c.push({name:i[h],data:o[i[h]]});t.chartConfig.series=c})},t.chartConfig={chart:{backgroundColor:null,zoomType:"x"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("ChartAreaCtrl",["$scope","$http","$rootScope",function(t,e,n){t.loadData=function(i){e.get("getDataFromFile?file="+i).success(function(e){e=e.data;var a=i.split("|");t.chartConfig.series=[];for(var o=[],s=0;s<a.length;s++){for(var l=[],r=0;r<e[a[s]].length;r++)"Date_AU_Global_FX"==i&&n.lastUpdate==e[a[s]][r].date&&(n.activeGlobalFXUsers=parseInt(e[a[s]][r][t.valueName])),"Date_AU_China_FX"==i&&n.lastUpdate==e[a[s]][r].date&&(n.activeChinaFXUsers=parseInt(e[a[s]][r][t.valueName])),"Date_AU_China_FC"==i&&n.lastUpdate==e[a[s]][r].date&&(n.activeChinaFCUsers=parseInt(e[a[s]][r][t.valueName])),"Date_AU_China_SC"==i&&n.lastUpdate==e[a[s]][r].date&&(n.activeChinaSCUsers=parseInt(e[a[s]][r][t.valueName])),"Date_NU_China"==i&&n.lastUpdate==e[a[s]][r].date&&(n.newChinaUsers=parseInt(e[a[s]][r][t.valueName])),"Date_NU_Global"==i&&n.lastUpdate==e[a[s]][r].date&&(n.newGlobalUsers=parseInt(e[a[s]][r][t.valueName])),l.push([Date.parse(e[a[s]][r].date),parseInt(e[a[s]][r][t.valueName])]);console.info(o),t.chartConfig.series.push({type:"area",name:t[a[s]],data:l})}console.info(t.chartConfig)})},t.chartConfig={chart:{backgroundColor:null,zoomType:"x"},xAxis:{type:"datetime",dateTimeLabelFormats:{month:"%b. %e",year:"%b"},title:{text:"Date"}},yAxis:{title:{text:"Users"},min:0},tooltip:{headerFormat:"<b>{series.name}</b><br>",pointFormat:"{point.x:%e. %b}: {point.y:.f} users"},legend:{enabled:!0},plotOptions:{area:{color:"#0038ff",marker:{enabled:null},lineWidth:1,states:{hover:{lineWidth:1}},threshold:null}},credits:{enabled:!1},series:t.chartSeries}}]),angular.module("inspinia").controller("BoxCtrl",["$scope","$http","$rootScope","$filter",function(t,e,n){n.lastUpdate=moment().subtract(1,"days").format("YYYY-MM-DD"),t.loadContestCount=function(){e.get("getContestCount").success(function(e){t.contentCount=e[0].count})}}]),$(document).ready(function(){function t(){var t=$("body > #wrapper").height()-61;$(".sidebard-panel").css("min-height",t+"px");var e=$("nav.navbar-default").height(),n=$("#page-wrapper").height();e>n&&$("#page-wrapper").css("min-height",e+"px"),n>e&&$("#page-wrapper").css("min-height",$(window).height()+"px"),$("body").hasClass("fixed-nav")&&$("#page-wrapper").css("min-height",$(window).height()-60+"px")}$(window).bind("load resize scroll",function(){$("body").hasClass("body-small")||t()}),setTimeout(function(){t()})}),$(function(){$(window).bind("load resize",function(){$(this).width()<769?$("body").addClass("body-small"):$("body").removeClass("body-small")})}),angular.module("inspinia").directive("sideNavigation",["$timeout",function(t){return{restrict:"A",link:function(e,n){e.$watch("authentication.user",function(){t(function(){n.metisMenu()})})}}}]).directive("minimalizaSidebar",["$timeout",function(t){return{restrict:"A",template:'<a class="navbar-minimalize minimalize-styl-2 btn btn-primary " href="" ng-init="minimalize()"><i class="fa fa-bars"></i></a>',controller:["$scope","$element",function(e){e.minimalize=function(){angular.element("body").toggleClass("mini-navbar"),!angular.element("body").hasClass("mini-navbar")||angular.element("body").hasClass("body-small")?(angular.element("#side-menu").hide(),t(function(){angular.element("#side-menu").fadeIn(500)},100)):angular.element("#side-menu").removeAttr("style")}}]}}]),angular.module("inspinia").run(["$templateCache",function(t){t.put("app/main/main.html",'<div class="wrapper wrapper-content animated fadeInRight"><div ng-include="\'components/box/sumbox.html\'"></div><div class="row chartPanel"><div ng-include="\'components/chart/chart_tu.html\'"></div></div><div class="row chartPanel"><div ng-include="\'components/map/map.html\'"></div><div ng-include="\'components/chart/chart_global_dist.html\'"></div></div><div class="row chartPanel"><div ng-include="\'components/chart/chart_fx_fc_sc.html\'"></div></div><div class="row chartPanel"><div ng-include="\'components/chart/overlap.html\'"></div><div ng-include="\'components/chart/veen.html\'"></div></div><div ng-controller="TabsDemoCtrl" class="tab-bar"><tabset><tab heading="Global"><div class="row chartPanel"><div ng-include="\'components/chart/chart_gb_nu.html\'"></div></div></tab><tab heading="China" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_nu.html\'"></div></div></tab></tabset></div><div ng-controller="TabsDemoCtrl" class="tab-bar"><tabset><tab heading="Global"><div class="row chartPanel"><div ng-include="\'components/chart/chart_gb_au.html\'"></div></div></tab><tab heading="China FX" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_au_fx.html\'"></div></div></tab><tab heading="China FC" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_au_fc.html\'"></div></div></tab><tab heading="China SC" ng-click="click()"><div class="row chartPanel"><div ng-include="\'components/chart/chart_cn_au_sc.html\'"></div></div></tab></tabset></div></div>'),t.put("components/box/sumbox.html",'<div class="row" ng-controller="BoxCtrl"><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-success pull-right">All Apps</span><h5>Total Users</h5></div><div class="ibox-content"><h1 class="no-margins" ng-bind="(($root.totalUsersChina + $root.totalUsersGlobal)|number)">0</h1><small>Users</small></div></div></div><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-info pull-right">Daily</span><h5>Active User</h5></div><div class="ibox-content"><h1 class="no-margins" ng-bind="(($root.activeChinaFXUsers + $root.activeChinaSCUsers + $root.activeChinaFCUsers + $root.activeGlobalFXUsers ) |number)">0</h1><div class="stat-percent font-bold text-success" ng-bind="(($root.activeChinaFXUsers + $root.activeChinaSCUsers + $root.activeChinaFCUsers + $root.activeGlobalFXUsers ) / $root.totalUsers * 100|number:2) + \'% active\'"><i class="fa fa-bolt"></i></div><small>Users</small></div></div></div><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-primary pull-right">Daily</span><h5>New Users</h5></div><div class="ibox-content"><h1 class="no-margins" ng-bind="(($root.newChinaUsers + $root.newGlobalUsers ) |number)">0</h1><div class="stat-percent font-bold text-navy" ng-bind="(($root.newChinaUsers + $root.newGlobalUsers ) / $root.totalUsers * 100|number:2) + \'% new\'"><i class="fa fa-bolt"></i></div><small>Users</small></div></div></div><div class="col-lg-3"><div class="ibox float-e-margins"><div class="ibox-title"><span class="label label-danger pull-right">Global</span><h5>Contest</h5></div><div class="ibox-content"><h1 class="no-margins" ng-init="loadContestCount()" ng-bind="(contentCount |number)">0</h1><small>Contests Held</small></div></div></div></div>'),t.put("components/chart/chart_cn_au_fc.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_China_FC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_China_FC=\'China Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_China_FC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_China_FC=\'China Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_China_FC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_China_FC=\'China Monthly Active Users\'"></div>'),t.put("components/chart/chart_cn_au_fx.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_China_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_China_FX=\'China Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_China_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_China_FX=\'China Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_China_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_China_FX=\'China Monthly Active Users\'"></div>'),t.put("components/chart/chart_cn_au_sc.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_China_SC=\'China Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_China_SC=\'China Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'China\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_China_SC=\'China Monthly Active Users\'"></div>'),t.put("components/chart/chart_cn_nu.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_NU_China\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily New Users\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Date_NU_China=\'China Daily New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_NU_China\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly New Users\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Week_NU_China=\'China Weekly New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_NU_China\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly New Users\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Month_NU_China=\'China Monthly New Users\'"></div>'),t.put("components/chart/chart_fx_fc_sc.html",'<div class="col-md-12" ng-controller="ChartCtrl"><highchart ng-init="loadData(\'Date_TU_China_FX|Date_TU_China_FC|Date_TU_China_SC\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'China FX FC SC\'"> <input type="hidden" ng-init="valueName=\'Total Users\'"> <input type="hidden" ng-init="Date_TU_China_FX=\'China FX\'"> <input type="hidden" ng-init="Date_TU_China_FC=\'China FC\'"> <input type="hidden" ng-init="Date_TU_China_SC=\'China SC\'"></div>'),t.put("components/chart/chart_gb_au.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_AU_Global_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Date_AU_Global_FX=\'Global Daily Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_AU_Global_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Week_AU_Global_FX=\'Global Weekly Active Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_AU_Global_FX\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly Active Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'Active Users\'"> <input type="hidden" ng-init="Month_AU_Global_FX=\'Global Monthly Active Users\'"></div>'),t.put("components/chart/chart_gb_nu.html",'<div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Date_NU_Global\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Daily New Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Date_NU_Global=\'Global Daily New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Week_NU_Global\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Weekly New Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Week_NU_Global=\'Global Weekly New Users\'"></div><div class="col-md-4" ng-controller="ChartAreaCtrl"><highchart ng-init="loadData(\'Month_NU_Global\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Monthly New Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Global\'"> <input type="hidden" ng-init="valueName=\'New Users\'"> <input type="hidden" ng-init="Month_NU_Global=\'Global Monthly New Users\'"></div>'),t.put("components/chart/chart_global_dist.html",'<div ng-controller="ChartComLineCtrl"><div class="col-md-6"><highchart ng-init="loadData(\'Date_Global_Country_TU|Global_Country_TU\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Global Distribution\'"></div></div>'),t.put("components/chart/chart_tu.html",'<div class="row" ng-controller="ChartCtrl"><div class="col-lg-12"><highchart ng-init="loadData(\'Date_TU_Global|Date_TU_China\')" config="chartConfig" class="span9"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'Total User Growth\'"> <input type="hidden" ng-init="valueName=\'Total Users\'"> <input type="hidden" ng-init="Date_TU_Global=\'Global\'"> <input type="hidden" ng-init="Date_TU_China=\'China\'"> <input type="hidden" ng-init="chartType=\'area\'"> <input type="hidden" ng-init="chartConfig.colors=[\'#0000ff\',\'ff0066\',\'2b908f\']"></div></div>'),t.put("components/chart/overlap.html",'<div class="col-md-6" ng-controller="ChartCtrl"><highchart ng-init="loadData(\'FX_FC_Overlap|FX_SC_Overlap|FC_SC_Overlap|FX_FC_SC_Overlap\')" config="chartConfig"></highchart><input type="hidden" ng-init="chartConfig.title.text=\'FX FC SC Users\'"> <input type="hidden" ng-init="chartConfig.subtitle.text=\'Overlaping Users\'"> <input type="hidden" ng-init="valueName=\'Common Users\'"> <input type="hidden" ng-init="FX_FC_Overlap=\'FX&FC\'"> <input type="hidden" ng-init="FX_SC_Overlap=\'FX&SC\'"> <input type="hidden" ng-init="FC_SC_Overlap=\'FC&SC\'"> <input type="hidden" ng-init="FX_FC_SC_Overlap=\'FX&FC&SC\'"></div>'),t.put("components/chart/veen.html",'<div class="col-md-6 venn-panel" ng-controller="VeenCtrl"><div id="venn" venn="vennData" ng-init="loadData(\'China_FX_FC_SC_TU\')"></div></div>'),t.put("components/common/content.html",'<div id="wrapper"><div ng-include="\'components/common/navigation.html\'"></div><div id="page-wrapper" class="gray-bg {{$state.current.name}}"><div ng-include="\'components/common/topnavbar.html\'"></div><div ui-view=""></div><div ng-include="\'components/common/footer.html\'"></div></div></div>'),t.put("components/common/footer.html",'<div class="footer"><div class="pull-right">FDT Dashboard v0.5</div><div><strong>Copyright</strong> FDT &copy; 2014-2015</div></div>'),t.put("components/common/ibox_tools.html",'<div class="ibox-tools dropdown" dropdown=""><a ng-click="showhide()"><i class="fa fa-chevron-up"></i></a> <a class="dropdown-toggle" href="" dropdown-toggle=""><i class="fa fa-wrench"></i></a><ul class="dropdown-menu dropdown-user"><li><a href="">Config option 1</a></li><li><a href="">Config option 2</a></li></ul><a ng-click="closebox()"><i class="fa fa-times"></i></a></div>'),t.put("components/common/navigation.html",'<nav class="navbar-default navbar-static-side" role="navigation"><div class="sidebar-collapse"><ul side-navigation="" class="nav metismenu" id="side-menu"><li class="nav-header"><div class="dropdown profile-element" dropdown=""><a class="dropdown-toggle" dropdown-toggle="" href=""><span class="clear"><span class="block m-t-xs"><strong class="font-bold">{{main.userName}}</strong></span> <span class="text-muted text-xs block">Example menu<b class="caret"></b></span></span></a><ul class="dropdown-menu animated fadeInRight m-t-xs"><li><a href="">Logout</a></li></ul></div><div class="logo-element">FDT</div></li><li ui-sref-active="active"><a ui-sref="index.main"><i class="fa fa-dashboard"></i> <span class="nav-label">User</span></a></li><li ui-sref-active="active"><a ui-sref="index.school"><i class="fa fa-graduation-cap"></i> <span class="nav-label">School</span></a></li><li ui-sref-active="active"><a ui-sref="index.region"><i class="fa fa-map-marker"></i> <span class="nav-label">Region</span></a></li></ul></div></nav>'),t.put("components/common/topnavbar.html",'<div class="row border-bottom"><nav class="navbar navbar-static-top white-bg" role="navigation" style="margin-bottom: 0"><div class="navbar-header"><span minimaliza-sidebar=""></span><form role="search" class="navbar-form-custom" method="post" action=""><div class="form-group"><input type="text" placeholder="" class="form-control" name="top-search" id="top-search"></div></form></div><ul class="nav navbar-top-links navbar-right"><li><span class="m-r-sm text-muted welcome-message" ng-bind="\'Last Updated Time : \' + $root.lastUpdate"></span></li><li><a href=""><i class="fa fa-sign-out"></i> Log out</a></li></ul></nav></div>'),t.put("components/map/map.html",'<div class="col-md-6 chartPanel" ng-controller="MapCtrl"><h3>Global Distribution on Map</h3><div id="map" map="" ng-init="data()"></div></div>'),t.put("components/school/contest.html",'<div class="col-lg-4" ng-controller="ContestCrol"><div class="ibox float-e-margins"><div class="ibox-title"><h5>School User In Contests</h5><div class="ibox-tools"><button type="button" class="btn btn-sm btn-primary">Export</button></div></div><div class="ibox-content"><div class="row"><div class="col-sm-3 m-b-xs product-option"><select class="input-sm form-control input-s-sm inline" ng-options="product.name for product in products" ng-model="selected" ng-change="changeProduct()"></select></div><div class="col-sm-3"></div></div><div class="table-responsive"><table class="table table-striped"><thead><tr><th class="hidden-xs">Contest Name</th><th class="hidden-xs"># Users</th></tr></thead><tbody><tr ng-repeat="contest in contestList"><td>{{contest.name}}</td><td>{{contest.count|number}}</td></tr></tbody></table></div></div></div></div>'),t.put("components/school/count.html",'<div class="row" ng-show="count_show"><div class="col-lg-3"><div class="widget style1"><div class="row"><div class="col-xs-4 text-center"><i class="fa fa-user fa-3x"></i></div><div class="col-xs-8 text-right"><span>Today Users</span><h2 class="font-bold" ng-bind="totalUser"></h2></div></div></div></div><div class="col-lg-3"><div class="widget style1"><div class="row"><div class="col-xs-4 text-center"><i class="fa fa-calendar-o fa-3x"></i></div><div class="col-xs-8 text-right"><span>Weekly Users</span><h2 class="font-bold" ng-bind="weekUser"></h2></div></div></div></div><div class="col-lg-3"><div class="widget style1"><div class="row"><div class="col-xs-4 text-center"><i class="fa fa-calendar fa-3x"></i></div><div class="col-xs-8 text-right"><span>Monthly Users</span><h2 class="font-bold" ng-bind="monthUser"></h2></div></div></div></div></div>'),t.put("components/school/school.html",'<div class="wrapper wrapper-content animated fadeInRight"><div class="row chartPanel" ng-controller="SchoolCtrl"><div class="row"><div class="col-md-3 selector pull-left"><ui-select ng-model="region.selected" on-select="onRegionSelect($item, $model)" theme="select2" class="form-control" title="Choose a Region"><ui-select-match placeholder="Select or search Region ...">{{$select.selected}}</ui-select-match><ui-select-choices repeat="item in regions | filter: $select.search"><div ng-bind-html="item | highlight: $select.search"></div></ui-select-choices></ui-select></div><div class="col-md-3 selector pull-left" ng-show="school_show"><ui-select ng-model="school.selected" on-select="onSchoolSelect($item, $model)" theme="select2" class="form-control" title="Choose a School"><ui-select-match placeholder="Select or search a school name ...">{{$select.selected.school_name}}</ui-select-match><ui-select-choices repeat="school in schools | filter: $select.search"><div ng-bind-html="school.school_name | highlight: $select.search"></div><small ng-bind-html="school.school_key | highlight: $select.search"></small></ui-select-choices></ui-select></div><div class="col-md-2 compare-switch" ng-show="school_show"><span style="float:right">Compare mode</span><div class="switch"><div class="onoffswitch"><input type="checkbox" ng-model="isCompard" class="onoffswitch-checkbox" id="compared"> <label class="onoffswitch-label" for="compared"><span class="onoffswitch-inner"></span></label></div></div></div></div><div ng-include="\'components/school/count.html\'"></div><div ng-include="\'components/school/user_growth.html\'"></div></div><div class="row chartPanel"><div ng-include="\'components/school/user_rank.html\'"></div><div ng-include="\'components/school/contest.html\'"></div></div></div>'),t.put("components/school/user_growth.html",'<div class="col-lg-12"><highchart config="chartConfig" class="span9"></highchart></div>'),t.put("components/school/user_rank.html",'<div class="col-lg-8" ng-controller="UserRankingCrol"><div class="ibox float-e-margins"><div class="ibox-title"><h5>School User Ranking(By ROI)</h5><div class="ibox-tools"><button type="button" class="btn btn-sm btn-primary">Export</button></div></div><div class="ibox-content"><div class="row"><div class="col-sm-3 m-b-xs product-option"><select class="input-sm form-control input-s-sm inline" ng-options="product.name for product in products" ng-model="selected" ng-change="changeProduct()"></select></div><div class="col-sm-3"></div></div><div class="table-responsive"><table class="table table-striped"><thead><tr><th></th><th>User Name</th><th class="hidden-xs">User Id</th><th class="hidden-xs">Register</th><th class="hidden-xs">Trade Count</th><th class="hidden-xs">Biggest Loss</th><th class="hidden-xs">Ranking</th><th>ROI</th></tr></thead><tbody><tr ng-repeat="userInfo in userRankList"><td><img class="img-circle" src="{{userInfo.serving_url || \'http://www.hkfdt.cn/favicon.ico\'}}"></td><td>{{userInfo.username}}</td><td>{{userInfo.user_id}}</td><td>{{userInfo.created|date}}</td><td>{{userInfo.tradescount|number}}</td><td class="neg-number">- {{userInfo.biggestloss|currency:"$"}}</td><td>{{userInfo.ranking}}</td><td>{{userInfo.unit_price|number:2}}</td></tr></tbody></table><div class="text-center"><div ng-hide="rankListLoaded"><i class="fa fa-spinner fa-spin"></i></div><button type="button" class="btn btn-sm btn-info" ng-class="userRankList.length>0?\'\':\'disabled\'" ng-click="loaded=true;loadMore()" ng-hide="loaded">Load More</button></div></div></div></div></div>')
}]);